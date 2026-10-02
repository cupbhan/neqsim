# 历史研究资料归档

**2026-10-02 后续接入：** 三个原件已派生出正常编译路径中的实现，并增加结果重算与物理验收接口。
范围和测试记录见[三项算法接入说明](../../../docs/development/shared-thermo-research-integration.md)。
下文保留归档时的说明；本目录的三个 `.java.txt` 原件仍不直接参与编译。

这里保存从 `neqsim-snapshot-backup-20260809` 提取的三份实验源代码和关键历史结果原件。原始路径、SHA256 和用途见上级目录的 [迁移清单](../migration-inventory.json) 与 [源码基线](../source-baseline.json)。文件保留原始字节；Git 属性禁止对这些原件转换换行符。

三个 Java 原件以 `.java.txt` 保存，尚未进入编译路径。模型 profile 用于按边界类型施加诊断性的水相互作用参数缩放；bifurcation corrector 是相分岔校正器；specified multiphase solver 用于指定相槽的两至四相闪蒸，本身不判断全局稳定性。它们需要各自的依赖核查和回归后才能成为公共内核能力。

`historical-` 文件是 2026 年 8 月已有研究的结果记录，供溯源和后续提取最小回归用例。本轮没有重新运行其完整 PVTsim 对标，也没有把这些汇总记录当成可独立重算的完整输入集。大型报告、图和原始扫描仍在原始备份目录，索引保留其路径和摘要。
