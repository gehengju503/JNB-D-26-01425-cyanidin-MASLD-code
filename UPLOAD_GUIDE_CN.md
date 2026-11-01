# JNB 修回代码上传详细步骤

## 一、上传前必须完成的核对

1. **不要直接公开当前文件夹中的个体级 Excel/CSV。** 代码仓库以脚本、变量字典、公开数据下载链接和可公开的汇总结果为主。若计划公开处理后的 NHANES 个体级数据，先核对 CDC/NCHS 数据使用与再分发要求，并由课题组负责人确认。
2. 在 R 4.3 或更高版本中运行 `R/00_install_packages.R`。把 `config.example.R` 复制为 `config.R`，填写本机输入路径。
3. 运行 `R/01_nhanes_analysis.R`，核对：总样本 3,882；CAP/LSM 3,875；HSI 3,647；FLI 1,756；Q4 对 Q1 的主模型与 BMI 敏感性模型结果是否逐项吻合正文和图注。
4. 核对纳排标准。当前可见代码尚未完整展示其他肝病、过量饮酒、妊娠、肝毒性药物等排除步骤。若论文声称执行了这些排除，必须补充可运行代码和每一步排除人数，再制作参与者流程图。
5. 核对 MASLD 算法。旧 Word 代码的代谢异常计数存在逻辑运算风险，不能直接上传。应由统计负责人确认腰围、血压、甘油三酯、HDL、糖代谢、HOMA-IR、CRP阈值及缺失值处理，并将最终算法写成单独函数和单元测试。
6. 运行 `R/02_go_kegg_enrichment.R`，确认输入基因表、背景基因集、BH校正、数据库版本和图表结果。已删除旧代码中随机生成 fold change 的做法，不能把随机值画成真实表达变化。
7. 按 `CODE_RESPONSE_MATRIX.csv` 补齐仍缺失的脚本。GSE202379/SRC/scTenifold主体代码已从RStudio缓存恢复并清理，但必须重新运行；仍需补齐原始质控/注释/donor-level差异分析、intervention RNA-seq、MR、PA/FFA CCK-8和saracatinib/TG分析。没有脚本或没有复跑时，不要在回复信中写“全部代码已公开”。
8. **必须重做Figure 9E-F和Figure 10。** 旧RStudio缓存含随机模拟火山图、人工GO表、错误修改下游靶基因及人为增加评分的代码。新包已删除这些步骤，旧图不能直接视为可复现结果。
9. 删除所有本机绝对路径、姓名、邮箱、口令、令牌、未公开受试者信息。确认 `config.R` 已加入 `.gitignore`。
10. 在一台不同电脑或新建 R 环境中从头复跑；保存 `sessionInfo.txt`，并让一名未参与写代码的合作者按 README 复现主要数字。

## 二、建立 GitHub 仓库

1. 登录 GitHub，点击右上角 `+` → `New repository`。
2. Repository name 建议：`JNB-D-26-01425-cyanidin-MASLD`。
3. Description 可写：`Reproducibility code for JNB-D-26-01425`。
4. 修回阶段可先选 `Private`；准备给编辑查看前，确认期刊是否接受匿名/私有链接。若需要公开 DOI，最终必须转为 `Public`。
5. 不勾选自动添加 README、LICENSE 或 `.gitignore`，因为本包已经包含这些文件。
6. 点击 `Create repository`。
7. 在本地打开本包根目录，依次执行 GitHub 页面给出的 “push an existing repository” 命令。典型命令为：

   ```text
   git init
   git add .
   git commit -m "Initial reproducibility release for JNB-D-26-01425"
   git branch -M main
   git remote add origin https://github.com/<组织或用户名>/JNB-D-26-01425-cyanidin-MASLD.git
   git push -u origin main
   ```

8. 在 GitHub 网页检查目录、README、CSV 是否正常显示；确认没有数据文件、`config.R`、缓存、临时图或账号信息。
9. 打开 `Settings` → `General` → `Social preview`（可选）与 `Features`；若需审稿匿名，不要在 README/提交作者中暴露不应披露的信息。
10. 补齐代码后打版本标签：`Releases` → `Draft a new release` → tag `v1.0.0` → title `JNB revision reproducibility release` → 发布。

## 三、用 Zenodo 生成永久 DOI

1. 登录 Zenodo，进入账户设置，在 `GitHub` 集成页连接 GitHub。
2. 找到刚建立的仓库，打开归档开关。
3. 回到 GitHub 发布 `v1.0.0` Release；Zenodo 会自动归档该版本并生成 DOI。
4. 在 Zenodo 记录页检查题名、作者、单位、ORCID、关键词、版本、许可证、funding 和 related identifiers。
5. 建议软件代码许可证用 MIT 或 BSD-3-Clause；数据许可证必须单独判断，不能用代码许可证覆盖受限数据。
6. 复制 Zenodo 的版本 DOI，而不是仅复制网页地址。把 DOI 写回 README 的 Suggested citation，然后再发布一个小版本（如 `v1.0.1`）以保存引用信息。

## 四、在稿件和回复信中更新

只有在相应脚本确已上传后，才使用下列措辞：

**Data availability 建议英文：**

> Public source data are available from NHANES, GEO (GSE48452 and GSE202379), eQTLGen, and FinnGen. Analysis code supporting the publicly reproducible components of this study is archived at Zenodo (DOI: [填写 DOI]) and mirrored at GitHub ([填写固定 release 链接]). Individual-level NHANES source data should be obtained from the official NHANES repository. The intervention RNA-seq data are [填写真实 accession；若尚未存储，不要虚构].

**Response to Reviewer 1, Comment 5 建议英文：**

> We have deposited the analysis code for [逐项列出实际已上传模块] in a versioned public repository and archived release (DOI: [填写 DOI]). The repository includes a README, software environment information, input specifications, and scripts implementing the reported filtering and multiple-testing procedures. The intervention RNA-seq data are available under accession [填写真实编号].

不要写笼统的 “all data and code are available”，除非矩阵中的每一项都已补齐并实际复跑。

## 五、期刊投稿系统上传

1. 登录 Elsevier Editorial Manager，打开修回任务，点击 `Revise Submission`。
2. 在 `Attach Files` / `Upload Files` 页面逐一上传：Response to Reviewers、带修订标记稿、清洁稿（若系统要求）、Supplementary Material、单独高清图、Highlights、Graphical Abstract，以及期刊要求的声明文件。
3. 代码通常不建议只作为一个无法更新的补充压缩包；优先填写 GitHub固定 Release 链接和 Zenodo DOI。若系统要求附件，可再上传本代码包 ZIP，并在说明中声明 DOI 为版本权威来源。
4. 文件类型要按系统选项准确指定；Response 不要误选为 Manuscript，Supplement 不要误选为 Figure。
5. 在 `View Submission` 生成合并 PDF 后逐页检查：公式、希腊字母、上下标、红色修订、图分辨率、补充表、链接和页码。
6. 在提交前最后核对正文 Data availability、回复信、README、Zenodo 记录中的 DOI/链接完全一致。
7. 提交后保存系统确认邮件、最终合并 PDF、GitHub release tag 和 Zenodo DOI 记录。
