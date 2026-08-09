# SAGD / NeqSim 迁移状态

更新日期：2026-08-08
分支：`sagd-integration`
上游基线：`875d1458f9346af78a68405ba4524370c2cbb66e`

## 目标

本目录是 NeqSim 的独立源码项目，用于维护可复用的物性、相平衡和 MCP 服务能力。SAGD 桌面软件只消费经过版本化、测试和校验的 JAR，不在产品仓库内维护 NeqSim Java 内核副本。

## 已完成

- 从 NeqSim 官方 Git 历史的精确提交建立干净仓库和独立分支。
- 从无 Git 元数据的旧快照迁入 SAGD 所需 runner、相包络求解代码及其最小生产依赖闭包。
- 排除旧快照中的 `target`、探针、报告和其他生成物。
- `mvnw.cmd spotless:check` 通过。
- `mvnw.cmd -DskipTests compile` 通过；3,100 个生产源文件成功编译。
- 将现场湿流体输入和 Ma et al. (2021) 拓扑基准资料固化为测试资源，消除对旧 `task_solve` 绝对路径的依赖。
- `FieldFluidRunnerTest`、`WaterIF97RunnerTest`、`HydrocarbonWaterBoundaryRegressionRunnerTest`、`HydrocarbonWaterModelIdentityTest` 和 `HydrocarbonWaterTopologyBenchmarkRunnerTest` 共 20 项全部通过。
- 扩展的 `HydrocarbonWater*Test` 数值/拓扑测试集 32/32 通过。
- 高含水状态点与两烃相包络已按不同组分基准计算：状态点保留全组分，包络使用归一化的非水组分，并在结果中披露排除的水相比例。
- MCP 组合构建使用未着色 NeqSim library 及其声明依赖，标准 NeqSim 发布仍保持 shaded 语义；重复 Java 类告警已经归零。
- 核心 Maven 工件固定输出时间并在所有打包步骤后规范化 ZIP；两个独立干净目录生成的普通 core library 均为 27,051,858 字节、SHA-256 `a3d8009d...66e2`，默认 shaded release 均为 57,519,244 字节、SHA-256 `f1a2e5fd...5975`，两种工件都直接字节一致。
- `neqsim-mcp-server` 固定 Quarkus 输出时间戳并在打包后规范化 ZIP 顺序与元数据；两个独立干净源码导出的 runner 均为 86,299,094 字节、SHA-256 `500c41af...2fc5`，直接字节比较一致。
- 该候选 JAR 已通过 SAGD 主项目真实 MCP/JSON 冒烟：TP flash、缓存、批处理、相包络质量门、黏度实验和 9 点 PVT 查表均通过。
- MCP 子项目已有 4 项直接契约测试，覆盖 TP flash、IF97、错误 JSON 边界和 `@Tool` 禁止 primitive 参数，4/4 通过。
- Java 版本边界已经统一：核心 `src/` 保持 Java 8 源码兼容，默认核心工件为 Java 17，MCP 是独立 Java 21 边界。
- 迁移候选的 98 个文件全部位于 Git index，工作区相对 index 无未暂存或未跟踪差异；`git diff --cached --check` 通过。
- 2026-08-08 使用固定 Adoptium 21.0.12+8 复验：Spotless 通过，核心目标测试 20/20 通过，MCP 契约测试 4/4 通过，并完成 core library 与 MCP runner 各自两次独立 clean build 同哈希和主项目真实客户端冒烟。

## 当前问题

当前工作树仍未形成迁移提交，因而可复现的 staged-tree JAR 仍不能冒充正式发布来源。服务字节级可复现和重复类治理已完成；剩余版本阻断是使用真实作者形成 clean commit，并从该提交生成、晋级带真实来源元数据的主项目运行 JAR。完整工具清单与传输层集成测试属于后续增强，不再阻断本轮构建基线。

## 下一步

1. 审核已暂存的 98 文件迁移集合，使用真实作者形成干净提交；不得虚构 Git 身份。
2. 从干净提交重跑双 clean build，生成版本化 JAR、SHA-256、许可证和来源清单，再由 SAGD 产品仓库晋级运行时并执行完整回归。
3. 扩展 MCP 完整工具清单与 HTTP/STDIO 传输层集成测试，并完成兼容性检查。
4. 完成 Visual Studio C++ Build Tools 安装后，关闭 SAGD Windows 原生模块构建门。

## 常用验收命令

```powershell
$env:JAVA_HOME = 'C:\mycodex\thermal wellbore simulator\.tools\jdk-21.0.12+8'
.\mvnw.cmd spotless:check
.\mvnw.cmd -DskipTests compile
.\mvnw.cmd '-Dtest=FieldFluidRunnerTest,WaterIF97RunnerTest,HydrocarbonWaterBoundaryRegressionRunnerTest,HydrocarbonWaterModelIdentityTest,HydrocarbonWaterTopologyBenchmarkRunnerTest' test
.\mvnw.cmd clean install -DskipTests -Dmaven.javadoc.skip=true -Dneqsim.shade.skip=true
.\mvnw.cmd -f neqsim-mcp-server\pom.xml -Plocal-dev spotless:check test
.\mvnw.cmd -f neqsim-mcp-server\pom.xml -Plocal-dev clean package -DskipTests -Dmaven.javadoc.skip=true
```
