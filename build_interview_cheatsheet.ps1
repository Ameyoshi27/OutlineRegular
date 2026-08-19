param(
    [string]$Source = 'E:\jt\面试\华测导航_三维重建算法工程师_一面冲刺复习.docx',
    [string]$Output = 'D:\outlineRegular\outlineRegular\华测导航_三维重建算法工程师_一面提词优化版.docx'
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
Add-Type -AssemblyName System.IO.Compression.FileSystem
Add-Type -AssemblyName System.IO.Compression

Copy-Item -LiteralPath $Source -Destination $Output -Force

$wordNs = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'
$xmlNs = 'http://www.w3.org/XML/1998/namespace'
$zip = [IO.Compression.ZipFile]::Open($Output, [IO.Compression.ZipArchiveMode]::Update)

function New-WordElement([xml]$Document, [string]$Name) {
    return $Document.CreateElement('w', $Name, $wordNs)
}

function Set-WordAttribute([xml]$Document, [Xml.XmlElement]$Element, [string]$Name, [string]$Value) {
    $attribute = $Document.CreateAttribute('w', $Name, $wordNs)
    $attribute.Value = $Value
    [void]$Element.Attributes.Append($attribute)
}

function Add-Run(
    [xml]$Document,
    [Xml.XmlElement]$Paragraph,
    [string]$Text,
    [bool]$Bold = $false,
    [string]$Color = '222222',
    [int]$Size = 21,
    [string]$Highlight = ''
) {
    $run = New-WordElement $Document 'r'
    $runProperties = New-WordElement $Document 'rPr'
    if ($Bold) {
        [void]$runProperties.AppendChild((New-WordElement $Document 'b'))
    }
    $colorNode = New-WordElement $Document 'color'
    Set-WordAttribute $Document $colorNode 'val' $Color
    [void]$runProperties.AppendChild($colorNode)
    $sizeNode = New-WordElement $Document 'sz'
    Set-WordAttribute $Document $sizeNode 'val' ([string]$Size)
    [void]$runProperties.AppendChild($sizeNode)
    $sizeCsNode = New-WordElement $Document 'szCs'
    Set-WordAttribute $Document $sizeCsNode 'val' ([string]$Size)
    [void]$runProperties.AppendChild($sizeCsNode)
    if ($Highlight) {
        $shading = New-WordElement $Document 'shd'
        Set-WordAttribute $Document $shading 'val' 'clear'
        Set-WordAttribute $Document $shading 'fill' $Highlight
        [void]$runProperties.AppendChild($shading)
    }
    [void]$run.AppendChild($runProperties)
    $textNode = New-WordElement $Document 't'
    $spaceAttribute = $Document.CreateAttribute('xml', 'space', $xmlNs)
    $spaceAttribute.Value = 'preserve'
    [void]$textNode.Attributes.Append($spaceAttribute)
    $textNode.InnerText = $Text
    [void]$run.AppendChild($textNode)
    [void]$Paragraph.AppendChild($run)
}

function New-Paragraph(
    [xml]$Document,
    [string]$Text,
    [ValidateSet('title','heading1','heading2','body','tip','danger','small')]
    [string]$Kind = 'body',
    [int]$Before = 0,
    [int]$After = 70,
    [bool]$KeepNext = $false
) {
    $paragraph = New-WordElement $Document 'p'
    $properties = New-WordElement $Document 'pPr'
    $spacing = New-WordElement $Document 'spacing'
    Set-WordAttribute $Document $spacing 'before' ([string]$Before)
    Set-WordAttribute $Document $spacing 'after' ([string]$After)
    Set-WordAttribute $Document $spacing 'line' '285'
    Set-WordAttribute $Document $spacing 'lineRule' 'auto'
    [void]$properties.AppendChild($spacing)
    if ($KeepNext) {
        [void]$properties.AppendChild((New-WordElement $Document 'keepNext'))
    }
    if ($Kind -eq 'tip' -or $Kind -eq 'danger') {
        $shading = New-WordElement $Document 'shd'
        Set-WordAttribute $Document $shading 'val' 'clear'
        Set-WordAttribute $Document $shading 'fill' $(if ($Kind -eq 'tip') { 'EAF3F8' } else { 'FFF2CC' })
        [void]$properties.AppendChild($shading)
        $indent = New-WordElement $Document 'ind'
        Set-WordAttribute $Document $indent 'left' '100'
        Set-WordAttribute $Document $indent 'right' '100'
        [void]$properties.AppendChild($indent)
    }
    [void]$paragraph.AppendChild($properties)

    switch ($Kind) {
        'title'    { Add-Run $Document $paragraph $Text $true '17365D' 34 }
        'heading1' { Add-Run $Document $paragraph $Text $true '17365D' 25 }
        'heading2' { Add-Run $Document $paragraph $Text $true '2F5597' 21 }
        'tip'      { Add-Run $Document $paragraph $Text $false '17365D' 20 }
        'danger'   { Add-Run $Document $paragraph $Text $false '7F6000' 20 }
        'small'    { Add-Run $Document $paragraph $Text $false '666666' 18 }
        default    { Add-Run $Document $paragraph $Text $false '222222' 20 }
    }
    return $paragraph
}

function New-LabelParagraph(
    [xml]$Document,
    [string]$Label,
    [string]$Text,
    [string]$Color = '17365D',
    [int]$After = 55
) {
    $paragraph = New-WordElement $Document 'p'
    $properties = New-WordElement $Document 'pPr'
    $spacing = New-WordElement $Document 'spacing'
    Set-WordAttribute $Document $spacing 'after' ([string]$After)
    Set-WordAttribute $Document $spacing 'line' '280'
    Set-WordAttribute $Document $spacing 'lineRule' 'auto'
    [void]$properties.AppendChild($spacing)
    [void]$paragraph.AppendChild($properties)
    Add-Run $Document $paragraph ($Label + '  ') $true $Color 20
    Add-Run $Document $paragraph $Text $false '222222' 20
    return $paragraph
}

function New-PageBreak([xml]$Document) {
    $paragraph = New-WordElement $Document 'p'
    $run = New-WordElement $Document 'r'
    $break = New-WordElement $Document 'br'
    Set-WordAttribute $Document $break 'type' 'page'
    [void]$run.AppendChild($break)
    [void]$paragraph.AppendChild($run)
    return $paragraph
}

try {
    $entry = $zip.GetEntry('word/document.xml')
    $reader = New-Object IO.StreamReader($entry.Open())
    try { [xml]$document = $reader.ReadToEnd() } finally { $reader.Dispose() }

    $namespaceManager = New-Object Xml.XmlNamespaceManager($document.NameTable)
    $namespaceManager.AddNamespace('w', $wordNs)
    $body = $document.SelectSingleNode('//w:body', $namespaceManager)
    $firstNode = $body.FirstChild
    $nodes = New-Object 'System.Collections.Generic.List[System.Xml.XmlNode]'
    function Add-Node([Xml.XmlNode]$Node) { [void]$nodes.Add($Node) }

    Add-Node (New-Paragraph $document '华测导航三维重建算法工程师｜一面屏幕提词版' 'title' 0 80 $true)
    Add-Node (New-Paragraph $document '使用方式：面试时主要停留在前 4 页；需要展开时，用 Word 导航窗格或 Ctrl+F 搜索关键词。后文保留完整复习资料。' 'small' 0 120)
    Add-Node (New-Paragraph $document '0. 先立住你的人设' 'heading1' 80 60 $true)
    Add-Node (New-Paragraph $document '测绘背景 + C++ 三维工程能力 + 点云几何处理与量产落地经验。当前强项不是从零研发 3DGS/SLAM，而是能把传感器数据、几何算法、软件平台和生产链路接起来。' 'tip' 0 90)
    Add-Node (New-LabelParagraph $document '岗位匹配' '强：点云处理、PCL/CGAL/OpenCV、C++/Qt/OpenGL、LOD 建模、ICP、R3LIVE 部署、多传感器标定、城市级生产。')
    Add-Node (New-LabelParagraph $document '主动补位' '3DGS、SfM/SLAM 后端优化、CUDA/Ceres/GTSAM/g2o 目前以原理理解为主；强调相关基础、学习路径和可迁移能力。')
    Add-Node (New-LabelParagraph $document '诚信边界' 'R3LIVE = 部署联调与标定，不说自研 SLAM；文物 = CloudCompare 粗配准/ICP/质检，不说自研 ICP；3DGS = 了解流程与关键问题，不说有项目经验。' '7F6000')
    Add-Node (New-Paragraph $document '90 秒自我介绍' 'heading2' 50 45 $true)
    Add-Node (New-Paragraph $document '面试官您好，我叫蒋涛，是西南交通大学测绘工程硕士，本科是地理信息科学，主要方向是点云软件开发和三维重建。我的核心经历是城市级建筑 LOD1.3 自动建模：我基于 C++、Qt、OpenGL 参与搭建点云处理与三维可视化平台，重点做建筑点云分层、轮廓提取、规则化建模和模块联调，项目已用于乐山约 2000 平方公里、1 万栋建筑生产。现在在天际航实习，处理 AI 建筑轮廓掩膜和倾斜摄影 OSGB 网格，打通采样、分栋裁剪、规则化建模和成果导出链路。我还做过 R3LIVE 手持设备部署、多传感器同步与外参标定，以及文物扫描数据的粗配准、ICP 精配准和质检。我和岗位最匹配的是点云几何处理、三维工程链路和生产落地；3DGS 与 CUDA 是我正在重点补齐的方向，希望能结合贵司手持、无人机和车载设备，把已有工程能力迁移到融合重建产品中。')
    Add-Node (New-Paragraph $document '回答节奏：结论一句 → 原理/做法两三句 → 项目证据一句 → 边界一句。每题控制在 40～90 秒，面试官追问再展开。' 'danger' 40 0)

    Add-Node (New-PageBreak $document)
    Add-Node (New-Paragraph $document '1. 项目高频短答｜先讲主项目' 'heading1' 0 65 $true)
    Add-Node (New-LabelParagraph $document 'LOD 流程' '机载点云预处理 → 建筑点提取/分栋 → 屋顶分层 → 二维投影与 AlphaShape 初始轮廓 → RANSAC 边线拟合 → 主方向/平行垂直约束规则化 → 层间拓扑处理与拉伸 → 闭合检查、交互修正和导出。')
    Add-Node (New-LabelParagraph $document '为什么规则化' '机载点云立面缺失、密度不均且边缘有噪声；直接曲面重建易鼓包、平滑棱角且面数高。生产需要规整、轻量、可编辑模型，因此采用数据驱动轮廓 + 建筑先验约束。')
    Add-Node (New-LabelParagraph $document 'AlphaShape' '能保留 L/U 形凹边界；α 太小会毛刺断裂，太大会抹平凹角，因此按点间距、建筑尺度生成多尺度候选并做质量筛选。')
    Add-Node (New-LabelParagraph $document 'RANSAC' '对离群点和局部缺失比最小二乘稳健；随机取最小样本拟合，按内点数评分，再用全部内点精化。项目中用于屋顶面/轮廓边线拟合。')
    Add-Node (New-LabelParagraph $document '评价指标' '几何：轮廓/平面/层高误差；正确性：漏建、错并、错分；拓扑：闭合、自交、重复面、非流形；生产：单栋耗时、人工修正率、批量稳定性。')
    Add-Node (New-LabelParagraph $document '真实难点模板' '先说具体失败现象，再定位是输入、参数、几何退化还是模块接口；展示中间结果；调整算法/阈值并用典型和失败样例回归。不要只说“调参数”。')
    Add-Node (New-Paragraph $document '三段次要经历｜每段 20～40 秒' 'heading2' 45 45 $true)
    Add-Node (New-LabelParagraph $document '天际航实习' '输入是 AI 建筑掩膜 + OSGB 网格；我参与网格采样、按轮廓分栋裁剪、规则化建模、成果导出及模块联调。重点是把上游识别结果转成稳定、可批处理的几何成果。')
    Add-Node (New-LabelParagraph $document 'R3LIVE' '做的是设备选型装配、Ubuntu/ROS 编译部署、时间同步、外参标定验证和彩色点云质检。理解紧耦合与传感器互补，但框架不是我研发。')
    Add-Node (New-LabelParagraph $document '文物配准' '人工粗配准让数据进入 ICP 收敛域，再做 ICP 精配准，最后看重叠区距离和边缘错位。我的贡献是实际数据处理和质量判断。')
    Add-Node (New-Paragraph $document '最可能被追问：你具体写了哪些模块/类？最难的一次问题？参数如何选？失败案例是什么？上线后指标如何验证？请只讲自己真正做过且能展开的内容。' 'danger' 40 0)

    Add-Node (New-PageBreak $document)
    Add-Node (New-Paragraph $document '2. JD 技术速查｜够答一面即可' 'heading1' 0 65 $true)
    Add-Node (New-LabelParagraph $document '3DGS 一句话' '用一组显式、可学习的各向异性三维高斯表示场景；将高斯投影到屏幕并按深度做 alpha 混合，训练快、渲染实时，但几何精度、大场景内存和稀疏视角仍是难点。')
    Add-Node (New-LabelParagraph $document '3DGS 流程' '多视图图像 → SfM 求相机位姿和稀疏点云 → 初始化高斯的位置/尺度/旋转/不透明度/颜色 → 可微 splatting 渲染 → 光度损失优化 → clone/split/prune 密度控制 → 实时新视角渲染。')
    Add-Node (New-LabelParagraph $document '为何需 SfM' '位姿决定每张图的投影关系，稀疏点云提供几何冷启动；位姿误差会直接造成重影、漂浮点和几何不一致。')
    Add-Node (New-LabelParagraph $document 'SfM / SLAM' '都估计相机运动和三维结构。SfM 偏离线高精度与全局 BA；SLAM 偏在线实时、关键帧、局部地图和回环，核心约束许多是共通的。')
    Add-Node (New-LabelParagraph $document 'SfM / MVS' 'SfM 先恢复相机位姿和稀疏结构；MVS 在位姿可靠的基础上恢复稠密深度/点云/表面。上游位姿误差会传到稠密重建。')
    Add-Node (New-LabelParagraph $document 'BA' '联合优化相机参数和三维点，以最小化所有观测的重投影误差；Ceres 常用于构建非线性最小二乘问题，工程上利用稀疏结构、局部 BA 或分块控制规模。')
    Add-Node (New-LabelParagraph $document 'ICP' '最近邻找对应 → SVD/最小二乘求刚体变换 → 更新 → 迭代收敛。依赖较好初值和足够重叠；低重叠、对称结构、噪声和动态物体会导致错误收敛。')
    Add-Node (New-LabelParagraph $document '融合关键点' '激光给尺度和几何，IMU 给高频运动与去畸变，相机给纹理和颜色。时间不同步会造成拖影/弯曲，外参不准会造成系统性错位；先查时间，再查坐标系和外参。')
    Add-Node (New-LabelParagraph $document '大场景性能' '先量化瓶颈，再做分块/瓦片、空间索引、减少复制、按需加载、视锥裁剪和 LOD；适合规则并行计算的热点再考虑 CUDA，不能把“上 GPU”当通用答案。')
    Add-Node (New-Paragraph $document '3DGS 细节追问时：先讲流程、表示、渲染和局限。具体迭代次数、阈值、学习率不是你的项目经验，不要抢答成固定标准；可以说“原论文常这样设置，工程中要随数据与实现验证”。' 'danger' 40 0)

    Add-Node (New-PageBreak $document)
    Add-Node (New-Paragraph $document '3. 临场话术｜卡住时也能稳住' 'heading1' 0 65 $true)
    Add-Node (New-LabelParagraph $document '不会某库' '这个库我还没有实际使用过。我理解它解决的是……，我做过的相近工作是……。如果进入项目，我会先跑通官方/项目最小示例，再用小数据验证输入输出与数值结果，最后接入现有链路并做性能分析。')
    Add-Node (New-LabelParagraph $document '不会 3DGS 实战' '我目前是原理和方案调研层面，没有把它包装成项目经验。我的优势是 SfM/点云/传感器/三维软件链路基础，以及处理工程数据和失败案例的经验；下一步会从 COLMAP + 原版 3DGS 跑通基线，再关注激光约束、几何正则和大场景分块。')
    Add-Node (New-LabelParagraph $document '问题不确定' '我先确认一下，您更关注算法原理，还是我在工程里会怎么实现和验证？')
    Add-Node (New-LabelParagraph $document '需要思考' '这个问题我先从输入、目标和约束拆一下……（停两秒）我的初步判断是……；其中我不确定的是……，会通过……验证。')
    Add-Node (New-LabelParagraph $document '纠正自己' '我刚才那句话不够准确，更严谨的表述应该是……。')
    Add-Node (New-Paragraph $document '为什么华测导航 / 为什么这个岗' 'heading2' 45 45 $true)
    Add-Node (New-Paragraph $document '这个岗位不是单一做论文指标，而是把手持、无人机、车载多传感器采集，接到 SfM/SLAM、点云优化、3DGS 表达和产品展示。我已有点云几何处理、手持融合系统部署、三维软件和城市级生产经验，和端到端工程链路很匹配；同时岗位里的 3DGS、优化后端和 GPU 也是我希望系统补齐的方向。')
    Add-Node (New-Paragraph $document '建议反问（选 2 个）' 'heading2' 45 45 $true)
    Add-Node (New-LabelParagraph $document '方向' '团队当前 3DGS 更偏视觉表达，还是也要求可量测的几何精度与测绘成果输出？')
    Add-Node (New-LabelParagraph $document '数据' '手持、无人机、车载数据目前如何融合？激光主要用于位姿约束、几何初始化，还是最终模型质量增强？')
    Add-Node (New-LabelParagraph $document '新人' '校招生进入后通常先负责哪一段：数据预处理/SfM-SLAM/稠密重建/3DGS 优化/工程加速？前三个月的评价标准是什么？')
    Add-Node (New-Paragraph $document '面试前 10 分钟：关掉消息弹窗；打开本文档与简历；测试摄像头/麦克风；纸上只写 4 个词：主项目、真实贡献、先结论、慢一点。' 'tip' 40 0)

    Add-Node (New-PageBreak $document)
    Add-Node (New-Paragraph $document '附录｜完整复习资料（需要展开时搜索）' 'heading1' 0 70 $true)

    foreach ($node in $nodes) {
        [void]$body.InsertBefore($document.ImportNode($node, $true), $firstNode)
    }

    $entry.Delete()
    $newEntry = $zip.CreateEntry('word/document.xml')
    $writerSettings = New-Object Xml.XmlWriterSettings
    $writerSettings.Encoding = New-Object Text.UTF8Encoding($false)
    $writerSettings.Indent = $false
    $stream = $newEntry.Open()
    $writer = [Xml.XmlWriter]::Create($stream, $writerSettings)
    try { $document.Save($writer) } finally { $writer.Dispose(); $stream.Dispose() }
}
finally {
    $zip.Dispose()
}

Write-Output $Output
