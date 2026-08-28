$ErrorActionPreference = 'Stop'

$docsDir = if ($PSScriptRoot) {
    $PSScriptRoot
} else {
    Join-Path (Get-Location) 'docs'
}
$docPath = Join-Path $docsDir 'maskonly_algorithm.docx'
$backupPath = Join-Path $docsDir 'maskonly_algorithm.before_update.docx'
$assetsDir = Join-Path $docsDir 'assets'
$overviewPath = Join-Path $assetsDir 'maskonly_overview.png'
$regularizationPath = Join-Path $assetsDir 'maskonly_regularization_flow.png'

New-Item -ItemType Directory -Path $assetsDir -Force | Out-Null
if (-not (Test-Path -LiteralPath $backupPath)) {
    Copy-Item -LiteralPath $docPath -Destination $backupPath
}

function Color([int]$r, [int]$g, [int]$b) {
    return $r + 256 * $g + 65536 * $b
}

function Add-FlowBox($slide, [double]$x, [double]$y, [double]$w, [double]$h,
                     [string]$title, [string]$body, [int]$fill, [int]$line,
                     [double]$titleSize = 17, [double]$bodySize = 11.5) {
    $body = $body.Replace('`r', [Environment]::NewLine)
    $shape = $slide.Shapes.AddShape(5, $x, $y, $w, $h)
    $shape.Fill.ForeColor.RGB = $fill
    $shape.Line.ForeColor.RGB = $line
    $shape.Line.Weight = 1.25
    $shape.Adjustments.Item(1) = 0.08
    $shape.TextFrame2.MarginLeft = 10
    $shape.TextFrame2.MarginRight = 10
    $shape.TextFrame2.MarginTop = 7
    $shape.TextFrame2.MarginBottom = 5
    $shape.TextFrame2.VerticalAnchor = 3
    $shape.TextFrame2.TextRange.Text = $title + "`r" + $body
    $shape.TextFrame2.TextRange.Font.NameFarEast = 'Microsoft YaHei'
    $shape.TextFrame2.TextRange.Font.Name = 'Arial'
    $shape.TextFrame2.TextRange.Font.Size = $bodySize
    $shape.TextFrame2.TextRange.Font.Fill.ForeColor.RGB = Color 31 41 55
    $shape.TextFrame2.TextRange.ParagraphFormat.Alignment = 2
    $shape.TextFrame2.TextRange.Paragraphs(1).Font.Bold = -1
    $shape.TextFrame2.TextRange.Paragraphs(1).Font.Size = $titleSize
    return $shape
}

function Add-Arrow($slide, [double]$x1, [double]$y1, [double]$x2, [double]$y2,
                   [int]$color = $(Color 71 85 105)) {
    $line = $slide.Shapes.AddConnector(1, $x1, $y1, $x2, $y2)
    $line.Line.ForeColor.RGB = $color
    $line.Line.Weight = 2
    $line.Line.EndArrowheadStyle = 3
    return $line
}

function Add-SlideTitle($slide, [string]$title, [string]$subtitle) {
    $titleBox = $slide.Shapes.AddTextbox(1, 36, 19, 888, 43)
    $titleBox.TextFrame2.TextRange.Text = $title
    $titleBox.TextFrame2.TextRange.Font.NameFarEast = 'Microsoft YaHei'
    $titleBox.TextFrame2.TextRange.Font.Size = 25
    $titleBox.TextFrame2.TextRange.Font.Bold = -1
    $titleBox.TextFrame2.TextRange.Font.Fill.ForeColor.RGB = Color 15 23 42
    $subBox = $slide.Shapes.AddTextbox(1, 38, 59, 884, 24)
    $subBox.TextFrame2.TextRange.Text = $subtitle
    $subBox.TextFrame2.TextRange.Font.NameFarEast = 'Microsoft YaHei'
    $subBox.TextFrame2.TextRange.Font.Size = 10.5
    $subBox.TextFrame2.TextRange.Font.Fill.ForeColor.RGB = Color 71 85 105
}

function Export-Flowcharts {
    $ppt = New-Object -ComObject PowerPoint.Application
    $ppt.Visible = -1
    $presentation = $ppt.Presentations.Add()
    $presentation.PageSetup.SlideWidth = 960
    $presentation.PageSetup.SlideHeight = 540

    try {
        $slide1 = $presentation.Slides.Add(1, 12)
        $slide1.Background.Fill.ForeColor.RGB = Color 248 250 252
        Add-SlideTitle $slide1 'Mask-only 总体处理流程' '从 GeoTIFF 掩膜到规则化 Shapefile；不依赖 OSGB 或外部点云'

        $blue = Color 219 234 254
        $green = Color 220 252 231
        $amber = Color 254 243 199
        $rose = Color 255 228 230
        $lineBlue = Color 37 99 235
        $lineGreen = Color 22 163 74
        $lineAmber = Color 217 119 6
        $lineRose = Color 225 29 72

        $b1 = Add-FlowBox $slide1 38 102 158 102 '1  掩膜分离' '近黑背景抑制`r颜色标签保留/二值回退`r腐蚀种子与连通域归属' $blue $lineBlue
        $b2 = Add-FlowBox $slide1 226 102 158 102 '2  栅格矢量化' '颜色漂移合并`r窄腰切分`rGDAL Polygonize' $blue $lineBlue
        $b3 = Add-FlowBox $slide1 414 102 158 102 '3  初始轮廓整理' '几何合并与包含清理`r保存 raw 初始轮廓`r拓扑保持平滑' $green $lineGreen
        $b4 = Add-FlowBox $slide1 602 102 158 102 '4  逐栋规则化' '统一方向检测`r拓扑主通道与定向重试`rVDP / StrictFallback 兜底' $amber $lineAmber
        $b5 = Add-FlowBox $slide1 790 102 132 102 '5  单栋验收' '有效性、面积、漂移`r墙面贴合、飞点`r方向合法性' $amber $lineAmber 16 10.5
        Add-Arrow $slide1 196 153 226 153 | Out-Null
        Add-Arrow $slide1 384 153 414 153 | Out-Null
        Add-Arrow $slide1 572 153 602 153 | Out-Null
        Add-Arrow $slide1 760 153 790 153 | Out-Null

        $b6 = Add-FlowBox $slide1 638 302 284 100 '6  建筑间重叠修复' '候选选择 → 带护栏的组平差 → 有界平移 → Difference 裁剪`r按规则化路径置信度和面积确定让位顺序' $rose $lineRose 17 11.5
        $b7 = Add-FlowBox $slide1 338 302 252 100 '7  写盘与重开审计' 'SyncToDisk 后关闭并重开图层`r再次修复/审计重叠`r删除面积 <15m² 或 bbox <20m² 的残片' $green $lineGreen 17 11.5
        $b8 = Add-FlowBox $slide1 38 302 252 100 '8  输出与调试产物' 'regularized_building.shp`r初始轮廓、方向系统、最佳假设`rraw 残差点等调试图层' $blue $lineBlue 17 11.5
        Add-Arrow $slide1 856 204 856 302 | Out-Null
        Add-Arrow $slide1 638 352 590 352 | Out-Null
        Add-Arrow $slide1 338 352 290 352 | Out-Null

        $note = $slide1.Shapes.AddTextbox(1, 38, 446, 884, 54)
        $note.TextFrame2.TextRange.Text = '关键边界：主路径失败后才进入有界重试；OBR 是最后的规则化降级输出，允许以 downgraded_quality 写出，并不等于通过完整质量门。'
        $note.TextFrame2.TextRange.Font.NameFarEast = 'Microsoft YaHei'
        $note.TextFrame2.TextRange.Font.Size = 12
        $note.TextFrame2.TextRange.Font.Fill.ForeColor.RGB = Color 127 29 29
        $note.Fill.ForeColor.RGB = Color 254 242 242
        $note.Line.ForeColor.RGB = Color 252 165 165

        $slide2 = $presentation.Slides.Add(2, 12)
        $slide2.Background.Fill.ForeColor.RGB = Color 248 250 252
        Add-SlideTitle $slide2 '单栋建筑规则化与降级路径' '一个基线 DirectionDetector 结果供各通道共享；失败重试可使用局部方向副本，不改写全局结果'

        $d = Add-FlowBox $slide2 40 102 196 86 'DirectionDetector' '边链化 → 固定带宽加权 KDE`r候选峰逐峰筛选 active`r单方向可做保守 PCA 纠偏' $blue $lineBlue 16 10.5
        $t = Add-FlowBox $slide2 286 102 214 86 '拓扑主通道' '格网吸附 → 共线合并`r转接构造 → 预检 → Ceres`r完整质量验收' $green $lineGreen 16 10.5
        $ok = Add-FlowBox $slide2 754 102 168 86 '通过' '保留固定拓扑`r进入建筑间重叠处理' $green $lineGreen 16 11
        Add-Arrow $slide2 236 145 286 145 | Out-Null
        Add-Arrow $slide2 500 145 754 145 | Out-Null

        $retry = Add-FlowBox $slide2 286 238 310 106 '原因驱动的有界重试' 'raw_kde：撤销不利的 PCA 纠偏`rstrict_geometry：加强合并/转接护栏`ralternate_direction：单方向失败时试强次峰或 PCA 轴' $amber $lineAmber 16 10.5
        Add-Arrow $slide2 393 188 393 238 | Out-Null
        Add-Arrow $slide2 596 291 704 291 | Out-Null
        $retryOk = Add-FlowBox $slide2 704 248 218 86 '任一重试通过' '使用该次局部方向上下文`r不重新运行全局方向检测' $green $lineGreen 15 10.5
        Add-Arrow $slide2 813 248 813 188 | Out-Null

        $vdp = Add-FlowBox $slide2 40 400 196 82 'VDP 备用简化' '生成低顶点数假设`r继承统一方向结果`r不直接作为最终输出' $blue $lineBlue 15 10.5
        $strict = Add-FlowBox $slide2 286 400 214 82 'StrictFallback' '优先多方向骨架`r再试单方向正交化' $amber $lineAmber 15 11
        $obr = Add-FlowBox $slide2 550 400 184 82 '方向 OBR' '最后规则形状`r可降级质量写出' $rose $lineRose 15 11
        $final = Add-FlowBox $slide2 784 400 138 82 '输出门槛' '简单有效多边形`r面积/bbox 下限' $green $lineGreen 15 10.5
        Add-Arrow $slide2 286 315 138 400 | Out-Null
        Add-Arrow $slide2 236 441 286 441 | Out-Null
        Add-Arrow $slide2 500 441 550 441 | Out-Null
        Add-Arrow $slide2 734 441 784 441 | Out-Null

        $slide1.Export($overviewPath, 'PNG', 2400, 1350)
        $slide2.Export($regularizationPath, 'PNG', 2400, 1350)
    }
    finally {
        try { $presentation.Close() } catch { }
        try { $ppt.Quit() } catch { }
        [System.Runtime.InteropServices.Marshal]::ReleaseComObject($presentation) | Out-Null
        [System.Runtime.InteropServices.Marshal]::ReleaseComObject($ppt) | Out-Null
    }
}

function Add-WordParagraph($doc, [string]$text, [int]$style = -1, [int]$alignment = 0) {
    $range = $doc.Range($doc.Content.End - 1, $doc.Content.End - 1)
    $range.Text = $text
    $range.Style = $style
    $range.ParagraphFormat.Alignment = $alignment
    $range.InsertParagraphAfter()
    return $range
}

function Add-WordImage($doc, [string]$path, [string]$caption) {
    $range = $doc.Range($doc.Content.End - 1, $doc.Content.End - 1)
    $image = $doc.InlineShapes.AddPicture($path, $false, $true, $range)
    $image.LockAspectRatio = -1
    $image.Width = 468
    $image.Range.ParagraphFormat.Alignment = 1
    $afterImage = $doc.Range($doc.Content.End - 1, $doc.Content.End - 1)
    $afterImage.InsertParagraphAfter()
    $captionRange = Add-WordParagraph $doc $caption -35 1
    $captionRange.Font.NameFarEast = 'Microsoft YaHei'
}

Export-Flowcharts

$word = New-Object -ComObject Word.Application
$word.Visible = $false
$word.DisplayAlerts = 0
$doc = $word.Documents.Open($docPath)

try {
    $doc.Content.Delete()
    $doc.PageSetup.PaperSize = 7
    $doc.PageSetup.TopMargin = 56.7
    $doc.PageSetup.BottomMargin = 56.7
    $doc.PageSetup.LeftMargin = 62.4
    $doc.PageSetup.RightMargin = 62.4

    $normal = $doc.Styles.Item(-1)
    $normal.Font.NameFarEast = 'Microsoft YaHei'
    $normal.Font.Name = 'Arial'
    $normal.Font.Size = 10.5
    $normal.ParagraphFormat.SpaceAfter = 6
    $normal.ParagraphFormat.LineSpacingRule = 1

    foreach ($styleId in @(-2, -3, -4)) {
        $style = $doc.Styles.Item($styleId)
        $style.Font.NameFarEast = 'Microsoft YaHei'
        $style.Font.Name = 'Arial'
        $style.Font.Color = Color 15 23 42
    }
    $doc.Styles.Item(-2).Font.Size = 16
    $doc.Styles.Item(-3).Font.Size = 13

    $title = Add-WordParagraph $doc 'Mask-only 模式算法说明' -63 1
    $title.Font.NameFarEast = 'Microsoft YaHei'
    $title.Font.Size = 24
    $title.Font.Color = Color 15 23 42
    $subtitle = Add-WordParagraph $doc 'outlineRegularTool：从 AI 掩膜到规则化建筑底面' -75 1
    $subtitle.Font.NameFarEast = 'Microsoft YaHei'
    $subtitle.Font.Color = Color 71 85 105
    Add-WordParagraph $doc '代码审核基线：2026-08-28 当前工作区。最近一次回归运行统计为 topology=210、VDP=3、StrictFallback=6、最终输出=218、自交=0；该统计仅用于回归对照，不是算法保证。' -1 1 | Out-Null

    Add-WordParagraph $doc '一、模式目标与边界' -2 | Out-Null
    Add-WordParagraph $doc 'Mask-only 模式只要求用户选择 AI 建筑掩膜 GeoTIFF 和输出 Shapefile，不读取 OSGB、点云或正射影像。程序从掩膜中提取建筑实例和初始轮廓，检测一个或多个方向系统，优先固定拓扑做联合平差，最后处理建筑间重叠。' | Out-Null
    Add-WordParagraph $doc '输入掩膜不要求“每栋建筑必定对应唯一颜色”。程序会抑制近黑背景；颜色数量可靠时保留颜色标签并合并轻微色漂，颜色过多时退化为二值前景；再利用腐蚀种子、连通域归属和窄腰切分恢复可能的建筑实例。' | Out-Null
    Add-WordParagraph $doc '结果目标是规则、有效且尽量贴合掩膜边界。拓扑通道及普通兜底结果必须通过相应质量检查；方向 OBR 是最后的规则化降级输出，代码允许以 downgraded_quality 写出，因此不能表述为“所有输出均通过完整质量门”。' | Out-Null

    Add-WordParagraph $doc '二、总体流程' -2 | Out-Null
    Add-WordImage $doc $overviewPath '图 1  Mask-only 模式总体处理流程'
    Add-WordImage $doc $regularizationPath '图 2  单栋建筑规则化、定向重试与降级路径'

    Add-WordParagraph $doc '三、初始轮廓提取' -2 | Out-Null
    Add-WordParagraph $doc '3.1 掩膜分离与矢量化' -3 | Out-Null
    Add-WordParagraph $doc '程序读取可用颜色波段和有效性掩膜，将近黑像素作为背景。颜色标签数量不超过可靠上限时保留实例颜色，并清理发丝状接缝、批次包裹环和轻微颜色漂移；超过上限时统一为二值前景。各颜色区域先腐蚀形成内部种子，再把原始前景像素归属到种子或原连通域，避免腐蚀直接缩小建筑。距离变换和分水岭用于切开窄腰粘连，最后由 GDALPolygonize 生成多边形。' | Out-Null
    Add-WordParagraph $doc '3.2 几何整理、raw 快照与平滑' -3 | Out-Null
    Add-WordParagraph $doc '矢量化后先做不依赖 OSGB 的几何合并、包含小块清理和窄颈切分。随后保存 initial_building_outline_raw.shp，作为未经平滑的像素边界证据。平滑阶段调用 OGR SimplifyPreserveTopology，容差按分辨率限制在约 0.20～0.45 m，并以有效性、面积变化、边界偏差和邻接关系作为护栏，生成 initial_building_outline.shp。' | Out-Null

    Add-WordParagraph $doc '四、统一方向检测' -2 | Out-Null
    Add-WordParagraph $doc '4.1 边链与加权 KDE' -3 | Out-Null
    Add-WordParagraph $doc 'DirectionDetector 在规则化前对平滑轮廓建立边链。近似同向的连续边被合并，链方向使用 TLS 主轴估计；弯曲度超限的链不被删除，而是降低长度权重，避免栅格楼梯边和渐变曲线主导结果。链方向折叠到 [0°,90°)，用固定 5° 带宽的圆环 KDE 累积密度。不同链通过长度和线性质量改变权重，当前代码没有“长边窄核、短边宽核”的自适应带宽。' | Out-Null
    Add-WordParagraph $doc '4.2 候选峰、生效方向与 PCA' -3 | Out-Null
    Add-WordParagraph $doc '局部密度极大值先经过 prominence 和最小角间距筛选，形成保留的候选方向系统；再按实际墙长、权重份额等证据逐峰确定 active。active 数量不少于 2 时 multiDirection=true。只有一个 active 系统且形状明显细长时，才允许在有限角差内用等弧长采样 PCA 主轴保守纠偏。debug_direction_systems.shp 保存的是经过候选筛选后保留的方向系统及 active/refined 状态，并不是所有原始 KDE 局部极大值。' | Out-Null
    Add-WordParagraph $doc '4.3 共享结果与局部重试' -3 | Out-Null
    Add-WordParagraph $doc '每栋建筑先生成一个基线 DetectedDirectionResult，拓扑通道、VDP 和 StrictFallback 共享该结果。首次拓扑成功时不会重新检测方向。失败后，raw_kde 和 alternate_direction 会复制方向上下文并在局部重试中覆盖主方向；它们不会回写或污染全局检测结果。因此“全程方向永不改变”不准确，正确语义是“基线只检测一次，失败重试可受控地使用局部方向副本”。' | Out-Null

    Add-WordParagraph $doc '五、拓扑保持规则化主通道' -2 | Out-Null
    Add-WordParagraph $doc '5.1 链吸附与候选拓扑' -3 | Out-Null
    Add-WordParagraph $doc '拓扑通道同时接收平滑初始轮廓、可关联的 raw 轮廓以及统一方向结果。边链被吸附到最近的 active 方向或其正交方向；相邻同轴链在整段偏差护栏允许时合并。相邻链优先用吸附线交点连接，异常交点则尝试方向合法的桥接、正交 dogleg、平行链 U-cap 或缺失封口恢复。转接构造不能引入方向系统外的自由斜边。' | Out-Null
    Add-WordParagraph $doc '5.2 双残差 Ceres 平差' -3 | Out-Null
    Add-WordParagraph $doc '候选通过自交、退化和面积等结构预检后进入 Ceres。优化变量主要是各吸附直线的偏移量，方向角保持锁定。raw 边界成功关联时，目标总权重为 raw 70%、smooth 30%；raw 点不足或关联不可靠时使用 smooth-only。旧文档把两者写反，且把 smooth 描述为主残差，均与当前代码不符。' | Out-Null
    Add-WordParagraph $doc '5.3 质量验收' -3 | Out-Null
    Add-WordParagraph $doc '验收覆盖简单多边形、面积比、质心漂移、墙面贴合、飞点、方向合法性和自交。预检只负责阻止明显结构错误，不能替代完整验收。固定拓扑候选失败时返回明确原因，供 wrapper 决定是否重试。' | Out-Null

    Add-WordParagraph $doc '六、失败重试与兜底' -2 | Out-Null
    Add-WordParagraph $doc '6.1 原因驱动的有界重试' -3 | Out-Null
    Add-WordParagraph $doc 'raw_kde 用于撤销造成自交等问题的 PCA 纠偏；strict_geometry 对合并覆盖区间、转接采样贴合和桥接几何启用更严格护栏；alternate_direction 在单方向候选持续失败且有足够证据时，尝试强次峰或高轴比 PCA 方向。重试次数有界，每次使用局部方向副本并重新完成构造和验收。' | Out-Null
    Add-WordParagraph $doc '6.2 VDP 与 StrictFallback' -3 | Out-Null
    Add-WordParagraph $doc '拓扑通道彻底失败后，VDP 生成较少顶点的备用假设并继承统一方向上下文。best hypothesis 是调试产物和 StrictFallback 的输入之一，不会作为未经规则化的最终轮廓直接写出。StrictFallback 在 active 系统不少于 2 时优先尝试多方向骨架，再尝试单方向正交化，最后生成主方向 OBR。Mask-only 最终路径不会直接回退并写出初始锯齿轮廓。' | Out-Null
    Add-WordParagraph $doc '6.3 OBR 的质量语义' -3 | Out-Null
    Add-WordParagraph $doc 'OBR 保证最后仍是方向规则的简单形状，但可能明显填平 U/L 形凹部或扩大面积。若完整质量检查不通过，当前代码仍可记录 [StrictOBR] accepted=downgraded_quality 后写出。因此应把它视为可定位、可统计的最后降级，而不是普通成功。' | Out-Null

    Add-WordParagraph $doc '七、建筑间重叠与最终审计' -2 | Out-Null
    Add-WordParagraph $doc '单栋结果写入前后均有重叠护栏。Mask-only 重叠修复按优先级选择候选，尝试带方向与位移约束的组 Ceres，再做有界平移；仍有重叠时对低优先级建筑执行 Difference，并检查有效性、面积和方向。写回后调用 SyncToDisk，关闭并重开输出图层，再运行一次最终重叠修复/审计。写出前和重开审计都会删除面积 <15 m² 或包围盒面积 <20 m² 的残片。' | Out-Null

    Add-WordParagraph $doc '八、局部曲线的当前状态' -2 | Out-Null
    Add-WordParagraph $doc '代码保留局部圆弧/椭圆弧检测与恢复接口，但 mask-only 当前设置 kMaskCurveDetectionDebugOnly=true，相关恢复调用不会执行。因此现阶段输出仍按直线方向系统规则化，不能把局部圆弧恢复描述为已启用能力；也没有启用贝塞尔曲线拟合。' | Out-Null

    Add-WordParagraph $doc '九、调试产物' -2 | Out-Null
    $tableRange = $doc.Range($doc.Content.End - 1, $doc.Content.End - 1)
    $table = $doc.Tables.Add($tableRange, 6, 2)
    $table.Borders.Enable = 1
    $table.Cell(1,1).Range.Text = '文件'
    $table.Cell(1,2).Range.Text = '当前含义'
    $table.Cell(2,1).Range.Text = 'initial_building_outline.shp'
    $table.Cell(2,2).Range.Text = '通过拓扑护栏平滑后的初始轮廓，供方向检测和规则化使用'
    $table.Cell(3,1).Range.Text = 'initial_building_outline_raw.shp'
    $table.Cell(3,2).Range.Text = '平滑前的像素级矢量轮廓，供 raw 残差关联与调试'
    $table.Cell(4,1).Range.Text = 'debug_direction_systems.shp'
    $table.Cell(4,2).Range.Text = '筛选后保留的候选方向系统；rank、active、raw_ang、refined 等字段描述是否生效及是否纠偏'
    $table.Cell(5,1).Range.Text = 'debug_best_hypothesis.shp'
    $table.Cell(5,2).Range.Text = 'VDP 最佳假设/兜底输入的调试几何，不代表一定被最终输出采用'
    $table.Cell(6,1).Range.Text = 'debug_mask_raw_residual_points.shp'
    $table.Cell(6,2).Range.Text = '成功关联并用于平差的 raw 残差采样点；未通过关联门的轮廓可能没有该图层记录'
    $table.Rows.Item(1).Range.Bold = -1
    $table.Rows.Item(1).Shading.BackgroundPatternColor = Color 219 234 254
    $table.Range.Font.NameFarEast = 'Microsoft YaHei'
    $table.Range.Font.Size = 9.5
    $table.Columns.Item(1).Width = 178
    $table.Columns.Item(2).Width = 318
    $endRange = $doc.Range($doc.Content.End - 1, $doc.Content.End - 1)
    $endRange.InsertParagraphAfter()

    Add-WordParagraph $doc '十、实现原则与已知限制' -2 | Out-Null
    Add-WordParagraph $doc '• 基线方向只检测一次；失败重试允许局部、可追踪的方向覆盖。' | Out-Null
    Add-WordParagraph $doc '• 候选峰不等于生效方向；multiDirection 由 active 系统数量派生。' | Out-Null
    Add-WordParagraph $doc '• 拓扑候选不允许方向系统外的自由边；构造无解时进入明确的降级路径。' | Out-Null
    Add-WordParagraph $doc '• raw 残差优先保存像素边界位置，但必须通过局部链关联，不能让无关边界点牵引整栋建筑。' | Out-Null
    Add-WordParagraph $doc '• OBR 降级、局部曲线暂停和少数 StrictFallback 案例是当前可见限制，应依靠路径统计和调试图层持续回归。' | Out-Null

    $section = $doc.Sections.Item(1)
    $header = $section.Headers.Item(1).Range
    $header.Text = 'outlineRegularTool  |  Mask-only 算法说明'
    $header.Font.NameFarEast = 'Microsoft YaHei'
    $header.Font.Size = 8.5
    $header.Font.Color = Color 100 116 139
    $header.ParagraphFormat.Alignment = 2
    $footer = $section.Footers.Item(1).Range
    $footer.Text = '第 '
    $footer.Collapse(0)
    $footer.Fields.Add($footer, 33) | Out-Null
    $footer.InsertAfter(' 页')
    $footer.Font.NameFarEast = 'Microsoft YaHei'
    $footer.Font.Size = 8.5
    $footer.Font.Color = Color 100 116 139
    $footer.ParagraphFormat.Alignment = 1

    $doc.Save()
}
finally {
    $doc.Close()
    $word.Quit()
    [System.Runtime.InteropServices.Marshal]::ReleaseComObject($doc) | Out-Null
    [System.Runtime.InteropServices.Marshal]::ReleaseComObject($word) | Out-Null
}

Write-Host "Updated: $docPath"
Write-Host "Flowcharts: $overviewPath; $regularizationPath"
