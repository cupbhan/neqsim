---
title: "Archived phase algorithms: compiled research APIs"
description: "Integration, acceptance checks, regression coverage and qualification limits for three archived phase algorithms."
---

# 三项历史研究算法接入

日期：2026-10-02。官方基线：v3.23.0。

原件保留在 [research-archive](../../distribution/cupbhan/research-archive/README.md)，不修改历史字节和 SHA256。
三个类已移入正常 Java 编译路径，必要数值回归纳入公共候选包和官方更新流程。
验证结果见[集成验收记录](../../distribution/cupbhan/research-integration/validation.json)。

## 公共接口

| 类 | 用途 | 入口 |
| --- | --- | --- |
| `HydrocarbonWaterBoundaryModelProfile` | 按烃—水边界族配置水的常数二元交互参数倍率 | 显式创建配置，通过 `createTemplate` 获得独立副本；默认倍率为 1 |
| `RetainedPhaseBifurcationCorrector` | 固定压力，以分岔方向的组成分离量为约束，联合校正温度、组成和新相比例 | `correct` 返回数值根和诊断 |
| `SpecifiedMultiphaseFlashSolver` | 固定温压，求解指定的两至四个相槽，允许重复相族 | `solve` 返回归一化组成、相比例和信赖域诊断 |

两项求解器提供共同的后续操作：

- `toThermodynamicSystem`：重建本实例产生的收敛结果，供诊断使用。
- `validateEquilibrium`：通过 `SpecifiedPhaseEquilibriumValidator` 独立检查物料平衡、归一化、逸度、相身份、相间距离和稳定性。
- `toValidatedThermodynamicSystem`：全部验收通过后才导出状态；失败抛出带原因的异常。

`isConverged` 只表示数值收敛。结果与产生它的求解器实例绑定，不能换用另一流体解释相组成数组。
三项研究接口不会自动接入默认闪蒸调用链；冻结的重油资源没有改变。

## 验收和模型资格

物料平衡和归一化残差须不超过 1e-8，对数逸度差不超过 1e-7；每个相的摩尔比例须大于 1e-10，
相间组成 L1 距离大于 1e-5。油/水身份通过 `LiquidPhaseClassification` 与个人 EOS 的质量贡献约定保持一致；
水摩尔分数高不自动意味着水相，气相身份由求值后的 EOS 状态决定。
适用的气、油、水 TPD 搜索须全部收敛，且未发现超出负稳定性阈值的非平凡分相趋势。

有限搜索不能证明数学意义上的全局稳定性，也不代替实验精度验证。
临界相合并、初生相、消失相交由专用相边界算法处理。四个相槽可求解不意味着任意流体存在稳定四相区。

参数配置仍标明 `experimentalCompatibilityProfile`、`requiresIndependentPhysicalValidation`，
并明确 `productionPromotionAllowed=false`。没有把历史 PVTsim 拟合倍率设为默认值，
也没有把软件间吻合当成实验标定。

历史水/轻质 TBP 常数交互参数调整也改为显式 `setLightPseudoWaterCompatibility(true)`，
普通流体沿用官方默认；该选项保留 CPA 参数、自定义重质 TBP 参数和原流体。
原因与回归见[合并后的行为审计](shared-thermo-compatibility-audit.md)。

## 回归与来源

覆盖副本隔离、独立边界族倍率、序列化、溢出拒绝；扰动初值后的两相/三相复算；
重复四相、消失相、错误相身份和额外分相趋势的拒绝；已知三相数学根与历史高压分岔反例；
结果归一化、流体绑定及调用者状态不变。旧高压四相扫描单独记录为较慢的研究回归。

新增源码和必要测试类列入[个人增强保护清单](../../distribution/cupbhan/personal-enhancements.json)，
统一构建必须确认这些测试实际执行。
[原迁移清单](../../distribution/cupbhan/migration-inventory.json)保留整理当时的状态，本页记录后续接入。
本轮不代表完整 PVTsim 重跑、实验标定或高压相边界网络闭合。

## rc.2 历史交付

研究接口首次交付时，本地公共运行版本切换为 `3.23.0-cupbhan.1-rc.2`，构建源码提交为
`f5db0c33b78d62da3a05130723430db44fff9595`。后续文档和运行指针提交不改变这个源码身份。
旧版本 `3.23.0-cupbhan.1-rc.1` 保留在版本目录中，可通过统一激活脚本回退。
候选包为本地可追溯发行包，未发布到公共 Maven 仓库。

| 验证范围 | 结果 |
| --- | --- |
| 核心选定回归，38 个测试类 | 190 通过，1 个已登记的上游禁用用例跳过，0 失败 |
| MCP 接口合同与真实 STDIO | 4 项测试、6 项运行检查通过，发现 76 个工具 |
| 三项研究算法专项 | 13 项通过，其中 11 项已包含在核心回归中，2 项历史扫描单独执行 |
| 仅加载发行 JAR 的 Python 调用 | 两相、三相、参数副本隔离、分岔数值根与相身份拒绝，4 项通过 |
| Python 绑定 | 13 项通过，覆盖 API 发现、显示、组分库和黏度 |
| thermal wellbore simulator | PVT、批处理、缓存、相包络验收、9 点查表、Python 调用及运行包身份检查通过 |
| T2WELL | 8 个 PR 状态点的 Fortran/NeqSim 对比通过；不同黏度关联式仅报告差异 |

验证范围是指定回归和实际调用检查，不等于全仓库测试、完整 UI 验收或实验标定。
本次未增加这三个研究接口的专用 MCP 工具或 GUI 控件；它们已可通过 Java 和 Python 调用。

修复了两项求解器在残差和雅可比计算时反复克隆 EOS、重复加载数据库的问题，
改为每次求解独享工作状态。相同分岔校正测试类本机用时从 640.0 秒降至 9.168 秒；
这是两次测试的观测值，不是统一平台性能基准。另有跨温压重复调用测试检查状态隔离。

完整版本身份、SHA256、命令和各平台验证见[交付记录](../../distribution/cupbhan/research-integration/delivery.json)，
候选构建原始记录见[候选包验收](../../distribution/cupbhan/research-integration/candidate-validation.json)。
构建记录中的 `productRuntimeSwitched=false` 是构建时状态；后续切换结果记录在交付记录中。

以后合并官方更新时，现有公共构建流程会检查个人增强源码和测试是否保留，并要求三项快速回归测试类实际执行。
重复相、消失相、错误油水身份和已发现额外分相趋势仍会阻止结果通过物理验收。
