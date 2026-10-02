"""Exercise the built shared engine over STDIO, with no third-party Python packages."""

import json
import queue
import subprocess
import threading
import time
from pathlib import Path


class McpClient:
    """A bounded line-delimited JSON-RPC client for the distribution smoke gate."""

    def __init__(self, jar, stderr_path, java="java"):
        self.messages = queue.Queue()
        self.sequence = 0
        self.stderr = open(stderr_path, "w", encoding="utf-8")
        self.process = subprocess.Popen(
            [java, "-Dquarkus.profile=stdio", "-jar", str(jar)],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=self.stderr,
            text=True, encoding="utf-8", errors="replace",
            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
        )
        self.reader = threading.Thread(target=self._read, daemon=True)
        self.reader.start()

    def _read(self):
        try:
            for line in self.process.stdout:
                try:
                    self.messages.put(json.loads(line))
                except json.JSONDecodeError:
                    continue
        finally:
            self.messages.put(None)

    def send(self, method, params=None, notification=False, timeout=180):
        message = {"jsonrpc": "2.0", "method": method}
        if params is not None:
            message["params"] = params
        if not notification:
            self.sequence += 1
            message["id"] = self.sequence
        self.process.stdin.write(json.dumps(message) + "\n")
        self.process.stdin.flush()
        if notification:
            return None
        deadline = time.monotonic() + timeout
        while True:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise TimeoutError("MCP response timeout: " + method)
            try:
                response = self.messages.get(timeout=remaining)
            except queue.Empty as exc:
                raise TimeoutError("MCP response timeout: " + method) from exc
            if response is None:
                raise RuntimeError("MCP process exited before responding: " + method)
            if response.get("id") != message["id"]:
                continue
            if "error" in response:
                raise RuntimeError(json.dumps(response["error"]))
            return response["result"]

    def call(self, name, arguments):
        result = self.send("tools/call", {"name": name, "arguments": arguments})
        if result.get("isError"):
            raise RuntimeError(json.dumps(result))
        for content in result.get("content", []):
            if content.get("type") == "text":
                data = json.loads(content["text"])
                if data.get("status") == "error" or data.get("error"):
                    raise RuntimeError(json.dumps(data))
                return data
        raise RuntimeError("No JSON text returned by " + name)

    def close(self):
        try:
            if self.process.poll() is None:
                self.process.stdin.close()
                try:
                    self.process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    self.process.terminate()
                    self.process.wait(timeout=10)
        finally:
            self.reader.join(timeout=2)
            self.process.stdout.close()
            self.stderr.close()


def run_smoke(jar, version, report_dir, repo_root, java="java"):
    """Validate actual runtime version, restored CPA selection and three-phase recovery."""
    report_dir = Path(report_dir)
    client = McpClient(jar, report_dir / "mcp-runtime.log", java)
    checks = []
    try:
        initialized = client.send("initialize", {
            "protocolVersion": "2025-11-25", "capabilities": {},
            "clientInfo": {"name": "cupbhan-shared-thermo-gate", "version": "1"},
        })
        actual_version = initialized["serverInfo"]["version"]
        if actual_version != version:
            raise RuntimeError(f"Runtime version mismatch: {actual_version} != {version}")
        checks.append({"name": "runtime-version", "passed": True, "version": actual_version})
        client.send("notifications/initialized", notification=True)
        tools = client.send("tools/list")
        names = {tool["name"] for tool in tools["tools"]}
        required = {"runFlash", "runFieldFluid", "runFluidFlash", "runWaterIF97"}
        policy_path = Path(repo_root) / "distribution/cupbhan/personal-enhancements.json"
        if policy_path.is_file():
            required.update(json.loads(policy_path.read_text(encoding="utf-8"))["requiredTools"])
        if not required.issubset(names):
            raise RuntimeError("Missing public tools: " + str(required - names))
        checks.append({"name": "tool-discovery", "passed": True, "count": len(names)})
        flash = client.call("runFlash", {
            "components": '{"methane":0.9,"ethane":0.1}',
            "temperature": 30, "temperatureUnit": "C", "pressure": 50,
            "pressureUnit": "bara", "eos": "PR", "flashType": "TP",
        })
        if not flash.get("fluid") or not flash.get("flash"):
            raise RuntimeError("PR flash omitted fluid or flash")
        checks.append({"name": "pr-flash", "passed": True})
        field = client.call("runFieldFluid", {"fieldFluidJson": json.dumps({
            "components": {"water": .9, "ammonia": .05, "CO2": .05},
            "temperatureC": 200, "pressureBara": 15, "temperatureMinC": 200,
            "temperatureMaxC": 200, "eos": "CPA", "reactive": False,
        })})
        if field.get("model") != "CPA" or field.get("speciation", {}).get("enabled") is not False:
            raise RuntimeError("CPA selected a different model or enabled reactive chemistry")
        checks.append({"name": "nonreactive-cpa-ammonia", "passed": True})
        payload_path = Path(repo_root) / "src/test/resources/neqsim/mcp/runners/field_wet_classic_payload.json"
        payload = json.loads(payload_path.read_text(encoding="utf-8"))
        payload.update(temperatureC=215, pressureBara=240, temperatureMinC=215,
                       temperatureMaxC=215, eos="SRK", reactive=False)
        wet = client.call("runFieldFluid", {"fieldFluidJson": json.dumps(payload)})
        if wet.get("flash", {}).get("numberOfPhases") != 3:
            raise RuntimeError("Water-rich recovery lost its three-phase solution")
        checks.append({"name": "water-rich-three-phase-recovery", "passed": True})
        water = client.call("runWaterIF97", {
            "temperatureK": 473.15, "pressureBara": 15, "temperatureMinK": 373.15,
            "temperatureMaxK": 473.15, "pointCount": 5,
        })
        if not water.get("fluid") or len(water.get("saturationCurve", [])) < 5:
            raise RuntimeError("IF97 response omitted state or curve")
        checks.append({"name": "water-if97", "passed": True})
        return {"passed": True, "checks": checks, "tools": sorted(names),
                "notCovered": ["HTTP transport", "full PVTsim rerun", "product UI",
                               "heavy-oil multimedia dedicated MCP endpoints"]}
    finally:
        client.close()
