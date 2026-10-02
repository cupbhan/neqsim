# 公共组分内核统一记录

本轮从现有增强主线整合 SAGD 独有改进，建立可重复执行的验证与候选包生成流程。记录日期为 2026 年 10 月 1 日，工作分支为 `integration/shared-thermo-20261001`。结果状态随实际检查更新。

## 原始基线

| 来源 | 固定版本 | 含义 |
| --- | --- | --- |
| 增强主仓库 | `e9b41932d6fa242022945394f47c340c65ae120d` | 本轮起点，POM 为 3.17.0，包含后续重油和多相修复 |
| 已合入的官方基线 | `66a866f1e2e82c2fbbd29a7af936e098cd29e849` | 本地主线已有的官方提交，不代表已查询到的最新官方版本 |
| SAGD 分支 | `760d0cb1933086c6c343b0ac57cf452e910513b9` | 基于 3.16.0，含独有 CPA 路由及氨组分修复 |
| 热井筒现用包 | SHA256 `83729788ec4f4dfa2d06026e006479b6c78c981be9ba60bd325d241d02bf5542` | 来源清单指向上述 SAGD 提交 |
| 无 Git 的研究快照 | `neqsim-snapshot-backup-20260809` | 研究资料和实验代码归档，不能分配一个虚构提交号 |

## 整合项目

| 项目 | 原始提交或来源 | 处理 |
| --- | --- | --- |
| SAGD 两相三相及 MCP 基础迁移 | `5d00dd3733` 与主线迁入记录 | 主线已有迁入版本，核对差异，避免整体覆盖 |
| 重油多介质模型与版本化资源 | `2d9b14dc4c` | 主线已有，纳入回归 |
| 高含水多相停滞恢复 | `a84fcf14c6` | 主线已有，纳入回归 |
| 普通 CPA 与电解质 CPA 区分、氨组分 | `760d0cb193` | 已通过保留原始提交来源的 cherry-pick 整合为 `570443c4d9` |
| 快照独有实验类 | 三个相边界实验类 | 已按原字节归档到 Git，保留摘要，尚未进入编译和公共接口 |
| 历史 PVTsim 对标资料 | 2026 年 8 月研究记录 | 建立带摘要的来源索引；本轮不把历史数值当作重新计算结果 |

## 需要保留的模型限制

高含水重油 200 至 300 bara 范围曾达到状态点求解回归通过，但完整边界网络仍未闭合。不能把收敛改善描述为该区间整体拓扑已经验证。

重油水 CO2 N2 NH3 的开发模型采用非反应 SRK、Peneloux 和温度相关经典混合参数。历史 NH3 剂量锚点为 0.1%、1% 和 5%，历史交付范围为 20 至 250 摄氏度、1 至 100 bara；不包含 pH、离子反应和反应热，也不代表现场标定。历史宽松交付标准通过，严格研究标准未通过，70 至 110 摄氏度气相消失边界仍有偏差。

同定义三相边界历史对标中，主边界平均绝对偏差约 0.878 摄氏度，次边界约 8.819 摄氏度。诊断性的水相互作用参数缩放未获生产晋级，不能作为默认模型参数。

## 执行状态

已创建整合分支，并从本机 SAGD 仓库读取原始提交到 `archive/sagd-integration` 跟踪引用。原主线与 SAGD 提交分别保存在 `archive/pre-unification-main-20261001` 和 `archive/pre-unification-sagd-20261001` 本地分支。核对最初 SAGD 迁移涉及的 94 个源码、测试及 MCP 文件，主仓库均有对应文件；这项检查表示迁移文件存在，不等于逐行或数值完全相同。

其中 87 个文件在忽略换行差异后内容一致。余下 7 个文件涉及后续闪蒸改进、氨组分别名、CPA 修复、MCP 工具与构建。原有 73 个 MCP public String 方法全部保留，新增 manageModel。细节见 [迁移清单](../../distribution/cupbhan/migration-inventory.json)。SAGD 工作副本的远端名已统一为 `origin` 指向个人 fork、`upstream` 指向官方；其源码提交和运行包未因此改变。

三份实验源码和六份关键历史结果原件已保存到 [研究归档](../../distribution/cupbhan/research-archive/README.md)，逐文件核对 SHA256，避免这些成果只存在于无 Git 的本机备份中。归档不代表实验算法已获工程验证。

2026 年 10 月 2 日首次运行核心 Spotless 和选定回归，75 项测试全部通过，无跳过。测试集合为 ComponentQueryTest、FieldFluidRunnerTest、WaterIF97RunnerTest、HeavyOilMultimediaFluidTest、TPmultiflashSolveStatusTest 和 HydrocarbonWater 系列测试。该结果不是全仓测试结论。

统一构建入口为 `devtools/build_shared_thermo.py`。它要求干净提交，使用个人版本分别构建 core library 和 MCP runner，检查内嵌版本，执行核心回归、MCP 契约测试和实际 STDIO 调用，最后生成来源、验证报告与校验清单。初始候选版本为 `3.17.0-cupbhan.1-rc.1`。脚本只生成本地候选目录；CI 工作流只保存构建附件，不发布或替换产品运行时。

发行流程的自动检查覆盖脏工作树、错误或旧的内嵌版本、缺失测试报告和失败或跳过的测试。首次全新候选构建中，核心通过，MCP 的独立 Spotless 检查发现 NeqSimTools 与 McpIdentityResolver 原有格式差异；已按该模块格式化并通过检查，排除注释和空白后的 Java token 序列未变化。相关 Java 兼容性技能也已纠正核心和 MCP 的版本边界，明确两个模块分别运行格式检查。

实际构建核对发现，Quarkus runner 不保留自身的 Maven pom.properties，其服务版本应从 JAR 清单的 Implementation-Version 读取，并通过 MCP initialize 再核对。内嵌核心仍要求 Maven 版本元数据，不能用服务版本代替。此差异已加入发行检查与测试。

候选目录已经生成，实际 runner 完成 6 项 STDIO 检查，并用热井筒原客户端再次验证闪蒸、缓存、批量闪蒸、相包络拒绝不合格结果、黏度实验和 9 点 CPA 查表。产品正式运行配置尚未切换。候选来源和完整验证记录以发行目录中的 SOURCE.json 和 validation.json 为准。

## 首个候选版本

版本为 `3.17.0-cupbhan.1-rc.1`，源码提交为 `bddee13a60be83aa72d157ede0b8cea23caae4b4`，候选标签为 `cupbhan-v3.17.0.1-rc.1`。本地发行目录为 `build/shared-thermo/releases/3.17.0-cupbhan.1-rc.1/`。测试结果及产物摘要另存为 [候选记录](../../distribution/cupbhan/candidates/3.17.0-cupbhan.1-rc.1.json)，便于在不提交大型 JAR 的情况下审查。

| 检查 | 本轮结果 |
| --- | --- |
| 核心指定回归 | 75 项通过，失败和跳过均为 0 |
| MCP 契约 | 4 项通过，失败和跳过均为 0 |
| 发行流程检查 | 10 项通过 |
| 实际 STDIO | 6 项通过，发现 74 个工具 |
| 热井筒原客户端 | 闪蒸、缓存、批量、质量拒绝、黏度实验及 9 点 CPA 查表通过 |
| 发行目录摘要 | 全部文件与 SHA256SUMS 一致 |

MCP runner 的 SHA256 为 `a1f58d255325294fe8c3513d4d7c431017ed4dd75312c888b45c759b1992ebfe`。这是候选版的局部回归结论；没有运行全仓测试、完整 PVTsim 重算或全部产品界面用例，也没有新增现场适用性声明。

重建需要 Python 3.11 或以上、Java 21、Git 和仓库自带的 Maven wrapper。在该标签的干净检出目录执行：

```text
python devtools/build_shared_thermo.py --version 3.17.0-cupbhan.1-rc.1
```

该命令已完整执行通过。脚本会拒绝覆盖已有同名发行目录；复现同一候选时使用新的干净检出目录。发行目录包含核心 JAR、依赖 POM、源码 JAR、MCP runner、LICENSE、源码基线、验证报告、SOURCE.json 和 SHA256SUMS。候选二进制目前保存在本机，没有发布到公共 Maven 或 GitHub Release。

核对接口时发现：Java 内核已经包含 HeavyOilMultimediaFluid，但 MCP 尚未提供产品客户端预留的 `runHeavyOilMultimediaFlash` 和 `runHeavyOilMultimediaBatch`。候选包必须披露这一缺口；不能把 Java 工厂测试通过当作这些专用接口可用。现有 `runFieldFluid` 与 `runFluidFlash` 的能力按其实际契约保留。

## 后续交付

首个候选内核完成后，热井筒、T2WELL、oilfield-calc-platform 分别锁定经过验证的版本。公共流体管理服务和 SDK 的独立仓库在内核基线稳定后建立。本轮不会把这些后续事项标记为已完成。
