library(shiny)
library(DT)
library(readr)
library(dplyr)
library(bregr)
library(visreg)
library(survival)
ui <- fluidPage(
  titlePanel("bregr 批量回归 · 全面 Shiny 模板"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "上传 CSV/TSV 数据", accept = c(".csv", ".tsv", ".txt")),
      checkboxInput("use_example", "没有数据？使用内置示例", value = TRUE),
      helpText("示例数据仅用于演示，系统会根据所选方法加载相应公开数据。"),
      hr(),
      uiOutput("method_ui"),
      uiOutput("y_block_ui"),
      uiOutput("x_ui"),
      uiOutput("x2_ui"),
      uiOutput("group_ui"),
      uiOutput("strata_ui"),
      uiOutput("xvar_ui"),
      numericInput("n_workers", "并行工作线程 n_workers", value = 1, min = 1, step = 1),
      actionButton("run", "开始建模", class = "btn-primary"),
      hr(),
      downloadButton("download_csv", "下载 Tidy 结果（CSV）")
    ),
    mainPanel(
      tabsetPanel(id = "tabs",
                  tabPanel("数据概览", DTOutput("head")),
                  tabPanel("Tidy 结果表", DTOutput("tbl")),
                  tabPanel("发表型森林图", plotOutput("forest", height = 520)),
                  tabPanel("环形森林图", plotOutput("forest_circle", height = 520)),
                  tabPanel("拟合线（visreg）", plotOutput("fitline", height = 520)),
                  tabPanel("残差诊断", plotOutput("residuals", height = 520)),
                  tabPanel("风险网络", plotOutput("risk_net", height = 620)),
                  tabPanel("列线图（Nomogram）", plotOutput("nomogram", height = 620)),
                  
                  # 这两个页签始终存在，但内容只在 method == 'coxph' 时显示
                  tabPanel("生存曲线（Cox）",
                           conditionalPanel("input.method == 'coxph'",
                                            plotOutput("surv", height = 520)
                           ),
                           conditionalPanel("input.method != 'coxph'",
                                            div(style = "padding:1rem;color:#666;",
                                                "提示：将方法切换为 CoxPH 才会绘制生存曲线。")
                           )
                  ),
                  tabPanel("Cox 诊断",
                           conditionalPanel("input.method == 'coxph'",
                                            plotOutput("cox_diag", height = 520)
                           ),
                           conditionalPanel("input.method != 'coxph'",
                                            div(style = "padding:1rem;color:#666;",
                                                "提示：将方法切换为 CoxPH 才会显示诊断图。")
                           )
                  )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # 1) 数据
  dataset <- reactive({
    if (!is.null(input$file)) {
      ext <- tolower(tools::file_ext(input$file$name))
      if (ext %in% c("csv","txt")) readr::read_csv(input$file$datapath, guess_max = 10000)
      else if (ext == "tsv") readr::read_tsv(input$file$datapath, guess_max = 10000)
      else validate("仅支持 .csv/.tsv/.txt")
    } else if (isTRUE(input$use_example)) {
      if (!is.null(input$method)) {
        if (input$method %in% c("coxph", "survreg", "cch")) {
          if (!requireNamespace("survival", quietly = TRUE))
            validate("请先安装 survival 包：install.packages('survival')")
          d <- survival::lung
          d$status <- ifelse(d$status == 2, 1L, 0L)  # 0=截尾,1=事件
          d$sex <- factor(d$sex, levels = c(1,2), labels = c("Male","Female"))
          d$ph.ecog <- factor(d$ph.ecog)
          d
        } else if (identical(input$method, "clogit")) {
          if (!requireNamespace("survival", quietly = TRUE))
            validate("请先安装 survival 包：install.packages('survival')")
          d <- survival::infert
          d$status <- as.integer(d$case)
          d$time <- 1
          d
        } else if (input$method %in% c("gamma", "inverse.gaussian")) {
          datasets::ChickWeight
        } else if (identical(input$method, "nls")) {
          datasets::Puromycin
        } else if (identical(input$method, "aov")) {
          datasets::PlantGrowth
        } else {
          mtcars
        }
      } else {
        mtcars
      }
    } else {
      validate("请上传数据或勾选示例数据")
    }
  })
  
  output$head <- renderDT({
    req(dataset())
    datatable(head(dataset(), 10), options = list(pageLength = 10), rownames = FALSE)
  })
  
  # 2) 方法
  output$method_ui <- renderUI({
    sel <- tryCatch(br_avail_methods(), error = function(e) c("gaussian","binomial","poisson","coxph"))
    method_desc <- c(
      gaussian = "连续因变量线性回归",
      binomial = "0/1 因变量逻辑回归",
      poisson  = "计数数据泊松回归",
      coxph    = "生存数据 Cox 比例风险模型",
      survreg  = "生存数据参数回归",
      clogit   = "匹配病例对照条件逻辑回归",
      cch      = "病例-队列设计生存模型"
    )
    desc <- method_desc[names(method_desc) %in% sel]
    help <- paste(names(desc), desc, sep = "：", collapse = "； ")
    tagList(
      selectInput("method", "建模方法", choices = sel, selected = sel[1]),
      helpText(help)
    )
  })
  
  # 3) 动态 UI
  observe({
    req(dataset())
    cols <- names(dataset())
    
    if (!is.null(input$method) && input$method %in% c("coxph","survreg","clogit","cch")) {
      output$y_block_ui <- renderUI({
        tagList(
          selectInput("y_time", "生存时间列（time）", choices = cols),
          selectInput("y_status", "结局状态列（0=截尾, 1=事件）", choices = cols)
        )
      })
      if (identical(input$method, "clogit")) {
        output$strata_ui <- renderUI({
          selectInput("y_strata", "匹配分层变量（strata）", choices = cols)
        })
      } else {
        output$strata_ui <- renderUI({ NULL })
      }
    } else {
      output$y_block_ui <- renderUI({
        selectInput("y", "因变量 Y", choices = cols)
      })
      output$strata_ui <- renderUI({ NULL })
    }
    
    output$x_ui <- renderUI({
      selectizeInput("x", "焦点自变量 X（批量）", choices = cols, multiple = TRUE)
    })
    output$x2_ui <- renderUI({
      selectizeInput("x2", "控制变量（可选）", choices = cols, multiple = TRUE)
    })
    output$group_ui <- renderUI({
      selectInput("group_by", "组变量（可选）", choices = c("无"="__none__", cols), selected = "__none__")
    })
    output$xvar_ui <- renderUI({
      selectInput("xvar", "拟合线横轴变量（visreg）", choices = cols)
    })
  })
  
  # 4) 建模
  breg_obj <- reactiveVal(NULL)
  
  observeEvent(input$run, {
    req(dataset(), input$method)
    dat <- dataset()
    group_by <- if (!is.null(input$group_by) && input$group_by != "__none__") input$group_by else NULL
    
    y_vec <- if (!is.null(input$method) && input$method %in% c("coxph","survreg","clogit","cch")) {
      req(input$y_time, input$y_status)
      st <- dat[[input$y_status]]
      if (!all(na.omit(unique(st)) %in% c(0,1)))
        showNotification("提醒：生存状态列应为 0/1（0=截尾,1=事件）。", type = "warning")
      if (identical(input$method, "clogit")) {
        req(input$y_strata)
        list(c(input$y_time, input$y_status), strata = input$y_strata)
      } else {
        c(input$y_time, input$y_status)
      }
    } else {
      req(input$y); input$y
    }
    
    validate(need(length(input$x) >= 1, "请至少选择一个自变量 X"))
    
    obj <- tryCatch({
      br_pipeline(
        data      = dat,
        y         = y_vec,
        x         = input$x,
        x2        = input$x2,
        method    = input$method,
        group_by  = group_by,
        n_workers = input$n_workers
      )
    }, error = function(e) {
      showNotification(paste("建模失败：", e$message), type = "error"); NULL
    })
    
    breg_obj(obj)
    if (!is.null(obj)) updateTabsetPanel(session, "tabs", selected = "Tidy 结果表")
  })
  
  # 5) 结果表 & 下载
  output$tbl <- renderDT({
    req(breg_obj())
    res <- br_get_results(breg_obj(), tidy = TRUE)
    datatable(res, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })
  
  output$download_csv <- downloadHandler(
    filename = function() paste0("bregr_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      req(breg_obj()); write.csv(br_get_results(breg_obj(), tidy = TRUE), file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )
  
  # 6) 图形
  output$forest <- renderPlot({ req(breg_obj()); br_show_forest(breg_obj()) })
  output$forest_circle <- renderPlot({ req(breg_obj()); br_show_forest_circle(breg_obj()) })
  
  output$fitline <- renderPlot({
    req(breg_obj(), input$xvar)
    if (!requireNamespace("visreg", quietly = TRUE)) {
      plot.new(); text(0.5, 0.5, "请安装 visreg 包以显示拟合线：install.packages('visreg')")
    } else {
      br_show_fitted_line(breg_obj(), xvar = input$xvar)
    }
  })
  
  output$residuals <- renderPlot({ req(breg_obj()); br_show_residuals(breg_obj()) })
  output$risk_net  <- renderPlot({ req(breg_obj()); br_show_risk_network(breg_obj()) })
  
  output$nomogram <- renderPlot({
    req(breg_obj())
    tryCatch({
      br_show_nomogram(breg_obj())
    }, error = function(e) {
      plot.new(); text(0.5, 0.5, paste("无法生成列线图：", e$message))
    })
  })
  
  output$surv <- renderPlot({
    req(breg_obj(), identical(input$method, "coxph"))
    br_show_survival_curves(breg_obj())
  })
  
  output$cox_diag <- renderPlot({
    req(breg_obj(), identical(input$method, "coxph"))
    br_show_coxph_diagnostics(breg_obj())
  })
}

shinyApp(ui, server)