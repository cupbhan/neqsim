# 本机 NeqSim 目录统一（2026-10-02）

唯一日常维护入口是 `C:\mycodex\neqsim`，分支为 `master`。此前的
`integration/shared-thermo-20261001` 已通过快进合入本地主线，不再另开一份源码维护。
官方跟踪远端仍为 `upstream`，个人增强仓库远端仍为 `origin`。

## 四个目录的比较与去向

| 原目录 | 实际角色 | 整理后的相对位置 |
| --- | --- | --- |
| `neqsim` | 3.17.0 基线及个人增强的 Git 主仓库 | 当前仓库根目录 |
| `neqsim-sagd` | 3.16.0 基线的旧产品分支，末次提交 `760d0cb193` | `.local-archive/20261002/neqsim-sagd/` |
| `neqsim-python` | Python/JPype 包装层，无独立 Git；原先附带旧 3.16.0 JAR | `bindings/python/` |
| `neqsim-snapshot-backup-20260809` | 无 Git 的研究快照，含源码、编译类、PVTsim 对比及研究产物 | `.local-archive/20261002/neqsim-snapshot-backup-20260809/` |

按文本统一换行后，SAGD 的 3,108 个生产 Java 文件中，2,978 个与当前主线相同，
130 个不同，没有主线缺失的生产 Java 文件。它的 1,631 个测试文件中 1,567 个相同、
64 个不同。差异结合提交历史判定，不能把较旧的不同文件覆盖到较新主线。
此前已迁入 SAGD 独有改进，具体见 `migration-inventory.json`。

快照有 3 个独有生产实验类、65 个独有测试文件，并有 129 个不同的生产 Java 文件。
完整快照保留原样；3 个实验类另已在 Git 中保留可审查的冻结副本，未直接启用为生产算法。
历史 `target/classes` 与旧 JAR 可能不一致，归档同时保存两者。

## 日常使用

```text
neqsim/
  src/                         Java 核心和数值测试
  neqsim-mcp-server/            MCP 接口
  bindings/python/             Python 包装层；共用同一增强内核
  distribution/cupbhan/         版本锁、来源、迁移和研究索引
  runtime/releases/<版本>/     本地固定运行包（忽略 Git）
  .local-archive/20261002/      历史目录和旧二进制（忽略 Git）
```

运行入口由 `distribution/cupbhan/runtime.json` 指定，当前为
`3.17.0-cupbhan.1-rc.1`。该包确实构建自 `bddee13a60be83aa72d157ede0b8cea23caae4b4`；
之后的目录、文档及包装层提交不会伪装成这个 JAR 的构建来源。
Python 包自身保留原接口版本 3.16.0，与所加载的增强内核版本分别管理。

安装 Python 接口（Java 21，Python 3.10+）：

```powershell
python -m venv bindings/python/.venv
bindings/python/.venv/Scripts/python.exe -m pip install -e ./bindings/python
```

仓库内 editable 安装可以自动找到根目录；其他位置的安装设置 `NEQSIM_HOME`
为本仓库或包含相同版本锁及运行包的部署目录。包装层验证 JAR 摘要及实际 Java 类来源，
避免已有 JVM 或外部 CLASSPATH 抢先加载另一个 NeqSim。

热井筒开发模式及其 `prepare:pvt` 使用同一版本锁；打包时复制经过校验的固定 JAR，
部署后不需要相邻源码目录。T2WELL 的 `tools/pvt/test_complib.py` 默认使用这个固定包，
仍保留 `--jar` 显式比较入口。oilfield-calc-platform 的比较脚本原本要求显式传入 JAR，
无需搬迁其源码。

## 后续升级

1. 只在此仓库的功能或同步分支中修改源码，按既定数值回归规则合入主线。
2. 从干净提交生成新版本候选，完成相应验证。
3. `python devtools/activate_shared_thermo.py build/shared-thermo/releases/<版本>`
   校验全部构建产物、安装不可覆盖的版本目录，并更新版本锁。
4. 验证使用平台，提交版本锁；热井筒执行 `npm run prepare:pvt` 更新部署副本。
   回退时激活保留的旧版本，不能覆盖旧版本同名 JAR。

## 保存与验证

搬迁前逐文件记录大小和 SHA-256，搬迁后对照校验。完整清单和结果保存在本地
`build/directory-consolidation/`，归档索引见同目录的 `move-manifest.json`。
三份搬迁目录合计 **57,912 个文件**，搬迁后全部与原始大小及 SHA-256 一致。
Python 包装层先以独立提交保存原始接口，再修改加载逻辑；其原始 JAR 另保留在
`.local-archive/20261002/python-original-runtime/`。
原 SAGD 的 `.git`、全部研究数据和 Python 原始捆绑 JAR 均保留。
大体积完整归档仅在本机，未声称已经上传 GitHub；关键来源和实验源码索引已受 Git 管理。

本次没有重跑完整 PVTsim 对标，也没有把实验模型提升为工程合格模型。
候选内核既有数值验证范围见 [统一记录](shared-thermo-unification.md)。

本次接入验证通过：6 项版本安装/摘要/不可覆盖规则测试、Python 原有 8 项接口测试、
实际 Python 闪蒸及 Java 类来源检查、旧 JVM 误加载拒绝检查、热井筒实际 MCP 闪蒸/批量/
缓存/黏度/9 点 PVT 表、适配器与固定包摘要检查，以及 T2WELL 的 PR 8 状态小规模对比。
T2WELL 这项检查用于验证新默认入口，不代替完整随机数值回归。
热井筒的 NeqSim 运行包来源一致性门已通过；整个软件发布的其他既有门槛不属于本次完成范围。
