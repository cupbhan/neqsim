---
title: "个人增强内核的官方版本同步"
description: "同步正式版、保留个人改进、运行回归并选择固定运行包的工作流"
---

# 官方更新与个人改进的持续合并

唯一维护仓库为 `cupbhan/neqsim`，主线为 `master`。同步目标默认是官方 GitHub
最新**正式发布版**，排除草稿和预发布；不将官方开发主线每个提交自动作为产品升级。
官方标签使用单独的 `refs/remotes/upstream-tags/` 跟踪，个人发行标签仍用 `cupbhan-v` 前缀。

## 日常流程

```text
官方正式版 → 检测与差异报告 → sync/upstream-v版本
                             ↓
                  三方合并、逐项解决冲突
                             ↓
             个人改进保留检查 + 数值与接口回归
                             ↓
                 固定候选包 + 来源与摘要
                             ↓
                    审阅并合入 master
                             ↓
                平台回归 → 激活运行包 → 可回退
```

源码同步和平台运行包选择是两个独立步骤。合并源码不会修改
`distribution/cupbhan/runtime.json`，也不会覆盖旧 JAR。

## 本地操作

在仓库根目录使用已选择的 Python 3.11+ 解释器运行：

```powershell
.\bindings\python\.venv\Scripts\python.exe devtools/sync_upstream.py
```

默认只查询正式版、取回其标签并生成差异报告。指定 `--tag v3.23.0` 可检查某个固定正式版。
报告位于 `build/upstream-sync/<标签>/report.json` 和 `REPORT.md`，包含官方提交、
个人起点、共同祖先、待合并提交数和双方均修改的文件。

确认工作树干净且当前位于待升级的主线后，增加 `--apply`：

```powershell
.\bindings\python\.venv\Scripts\python.exe devtools/sync_upstream.py --apply
```

脚本保存旧提交引用并创建 `sync/upstream-v版本`，使用 Git 三方 merge 保留双方历史。
遇到冲突即停在可检查的状态，不猜测 ours/theirs，不重置本地修改，不覆盖已有同步分支。
再次执行同一版本的合并会报告已集成；已有但未完成的同步分支需要继续处理该分支。
解决冲突后分别执行 core/MCP 的 Spotless、提交，并按报告的正式版 SHA 更新
`distribution/cupbhan/upstream.json`。选择 `git merge --abort` 可撤回尚未提交的合并；
已经提交的合并通过普通 revert 回退，不能强推覆盖主线。

## GitHub 工作流

入口是 Actions 中的 **Check and prepare official NeqSim updates**：

- 每周一 07:23 UTC 自动检测，生成差异报告与运行摘要。
- 手动 **Run workflow**，保留 `prepare=true`，可以合并新分支并构建验证候选包；
  `tag` 留空即使用最新正式版。
- 只有验证通过才推送同步分支；运行摘要给出比较和创建 PR 的链接。
- 冲突、测试失败或版本来源不一致会保留诊断附件；不会推送未验证的同步分支到主线。
- 当前仓库不允许 Actions 机器人创建/批准 PR，因此流程使用比较链接，未扩大仓库权限。
- 候选构建直接在此工作流内执行，避免依赖 `GITHUB_TOKEN` 推送触发另一条工作流。

`shared_thermo_candidate.yml` 仍支持 PR 和手动候选验证。版本从当前官方基线生成，
不再固定为 3.17.0；CI 候选号带运行 ID。CI 上传附件，不自动发布 Maven 或正式 Release。
官方继承的发布工作流仍限制在 `equinor/neqsim` 运行。

## 个人改进的保留与验证

`distribution/cupbhan/personal-enhancements.json` 是每次同步的兼容性清单：

- 保留重油多介质模型、含水相边界研究模块及关键接口。
- 冻结三份重油参数资源的 SHA-256；更新参数必须明确审查，不能随合并静默漂移。
- 保留原有 74 个 MCP 工具名称；构建后的实际工具发现也检查该集合。
- 强制执行原有 19 个核心测试套件，同时覆盖新版官方的水富集、反应状态和相包络回归。
- 拒绝缺失测试报告、失败、未登记的跳过和版本来源不匹配。
  官方原有的纯甲烷相包络禁用用例单独登记并如实报告，不计为通过；新增跳过会阻止候选发布。

常规候选门槛采用清单中的具体相包络测试套件。`PTPhaseEnvelopeTest` 中的 UMR-PRU
大规模温压粗网格扫描属于额外集成验证，本次扩展尝试在该扫描期间中止，未记为通过；
Michelsen、连续性、参考 EOS、截断和稳健性检查仍在常规门槛中。

个人改进已由官方实现替代时，记录替代原因，并用原回归及官方回归共同验证。
不要求旧实现逐字不变，也不能仅凭 Git 自动合并成功判定数值兼容。

`upstream.json` 记录**当前**官方基线。`source-baseline.json` 保持原历史基线及
PVTsim 记录原义；候选 `SOURCE.json` 同时记录当前官方提交与实际构建提交。

## 发布与回退

使用 `devtools/build_shared_thermo.py --version <官方版本>-cupbhan.<序号>-rc.<序号>`
从干净提交构建，版本不得复用。验证报告、源码 JAR、核心 JAR、MCP runner、许可证、
当前官方来源、个人改进清单和 SHA-256 一起保存。

完成平台回归后，用 `devtools/activate_shared_thermo.py` 指向该候选目录，提交版本锁。
热井筒执行 `npm run prepare:pvt` 和相应验证，更新其部署副本及来源报告。
Python 和 T2WELL 比较入口读取同一版本锁。回退时激活保留的旧发行目录，
再重新准备平台部署副本，不修改旧版本文件。

完整 PVTsim 对标、现场验证及全部产品 UI 验收是独立门槛，不能由本流程的局部回归代替。

## 本次同步

本次目标是官方 `v3.23.0`，提交 `280152d97763672fb79c9139e421adf85663a930`，
个人起点为 `3e9d6b4214eceeba26d77ff3808997ff19ce347b`。官方新增 1,002 个提交，
双方同时修改 20 个文件，实际产生 7 个文件冲突。执行结果记录在
`distribution/cupbhan/upstream-sync/`，原始诊断日志在 `build/upstream-sync/`。

相关说明：[仓库职责](shared-thermo-repository.md)、[目录统一](local-directory-consolidation.md)、
[官方 3.23.0 发布说明](https://github.com/equinor/neqsim/releases/tag/v3.23.0)、
[GitHub 工作流触发规则](https://docs.github.com/en/actions/how-tos/write-workflows/choose-when-workflows-run/trigger-a-workflow)。
