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
默认闪蒸调用链和冻结的重油资源没有改变。

## 验收和模型资格

物料平衡和归一化残差须不超过 1e-8，对数逸度差不超过 1e-7；每个相的摩尔比例须大于 1e-10，
相间组成 L1 距离大于 1e-5。油/水身份沿用烃—水研究路径的水摩尔分数 0.5 分类约定。
适用的气、油、水 TPD 搜索须全部收敛，且未发现超出负稳定性阈值的非平凡分相趋势。

有限搜索不能证明数学意义上的全局稳定性，也不代替实验精度验证。
临界相合并、初生相、消失相交由专用相边界算法处理。四个相槽可求解不意味着任意流体存在稳定四相区。

参数配置仍标明 `experimentalCompatibilityProfile`、`requiresIndependentPhysicalValidation`，
并明确 `productionPromotionAllowed=false`。没有把历史 PVTsim 拟合倍率设为默认值，
也没有把软件间吻合当成实验标定。

## 回归与来源

覆盖副本隔离、独立边界族倍率、序列化、溢出拒绝；扰动初值后的两相/三相复算；
重复四相、消失相、错误相身份和额外分相趋势的拒绝；已知三相数学根与历史高压分岔反例；
结果归一化、流体绑定及调用者状态不变。旧高压四相扫描单独记录为较慢的研究回归。

新增源码和必要测试类列入[个人增强保护清单](../../distribution/cupbhan/personal-enhancements.json)，
统一构建必须确认这些测试实际执行。
[原迁移清单](../../distribution/cupbhan/migration-inventory.json)保留整理当时的状态，本页记录后续接入。
本轮不代表完整 PVTsim 重跑、实验标定或高压相边界网络闭合。
