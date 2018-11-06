#CURRENT APPS:
#LPA - LPA Tribolium model
#LVcomp - Lotka-Volterra competion model
#LVpred - Lotka-Volterra predation model
#NB - Nicholson-Bailey Host-Parasitoid model
#May - May et al's Negative-Binomial parasitoid-host model
#Ricker - The Ricker ("discrete logisitc model")
#RM - Rosenzweig-MacArthur Predator-Prey model
#SEIR - the seasonally forced SEIR model
#SEIRS - the unforced SEIRS model
#SIR - the unforced SIR model
#TSIR - the unforced TSIR model with demographic an environmental stochasticity
#
#TO RUN: 
#source("allapps3.R")
#runApp(APP)

#require(shiny)
#require(scatterplot3d)
#require(deSolve)
#require(phaseR)


#' display a number
#' @param x A number
#' @return y A number of class "test"
#' @examples
#' test(4)
#' @export
test<-function(x){y=x
class(y)="test"
 return(y)}
 
 #' print a number of class "test"
 # @param y an object of class "test"
 #' @export
 print.test=function(y){
 cat("The value is ", y)
 }

Attr=c("ATTRIBUTION: This App was written by Ottar N. Bjornstad (onb1@psu.edu) and is licensed under the Creative Commons attribution-noncommercial license (http://creativecommons.org/licenses/by-nc/3.0/). Please share & remix non-commercially, mentioning its origin.")

#' Launch a shiny-app simulating the LPA model
#' @details
#' Launch app for details
#' @examples
#' \dontrun{LPA.app}
#' @export
LPA.app=shinyApp(
# This creates the LPA User Interface (UI)
ui = pageWithSidebar(
  tags$head(tags$style(
    HTML('
         #sidebar1 {
            background-color: #ECECEC;
        }
    
    #sidebar2 {
      background-color: #ECECEC
    }')
  )),
titlePanel("LPA model"),
fluidRow(
column(4, id = "sidebar2",
fluidRow(column(5, id = "sidebar1",
sliderInput("b", "b:", 6.598,
              min = 0, max = 10, step=0.1),
sliderInput("cel", "cel:", 1.209e-2,
              min = 0, max = 1, step=0.01),
sliderInput("cea", "cea:", 1.155e-2,
              min = 0, max = 1, step=0.01)
),
column(5, offset = 1, id = "sidebar1",
sliderInput("cpa", "cpa:", 4.7e-3,
              min = 0, max = 1, step=0.01),
sliderInput("mua", "mua:", 7.629e-3,
              min = 0, max = 1, step=0.01),
sliderInput("mul", "mul:", 0.2055,
              min = 0, max = 1, step=0.01)
#numericInput("L0", "initial L:", 50,
#              min = 0, max = 100),
#numericInput("P0", "initial P:", 0,
#              min = 0, max = 100),
#numericInput("A0", "initial A:", 0,
#              min = 0, max = 100)
),
column(1)),
fluidRow(
column(6, offset = 3, id = "sidebar1",
sliderInput("T", "Time range:",
                  min = 1, max = 500, value = c(1,500)),
checkboxInput("li", "lines", TRUE)),
column(3))
),
column(8, 
mainPanel(
  tabsetPanel(
      tabPanel("Time", plotOutput("plot1")), 
      tabPanel("2D Phase plane", plotOutput("plot2")),   
      tabPanel("3D Phase plane", plotOutput("plot3")),
      tabPanel("Details", 
           withMathJax(
            helpText("MODEL:"),
            helpText("Larvae $$L_{t+1} = b A_t e^{-c_{ea}A_t}e^{-c_{el}L_t}$$"),
           helpText("Pupae $$P_{t+1} = (1-m_l)L_t$$"),
           helpText("Adults $$A_{t+1} = P_te^{-c_{pa}A_t} + (1-m_a)A_t$$"),
           helpText("REFERENCE: Costantino RF, Desharnais RA, Cushing JM, Dennis B (1997) Chaotic dynamics 
            in an insect population. Science 275: 389-391"),
           helpText(eval(Attr))
           )))
)))
),

# This creates the 'behind the scenes' code (Server)
server <- function(input, output) {

lpa = function(b=6.598, cel=1.209e-2, cea=1.155e-2, mul=0.2055, 
      cpa=4.7e-3, mua=7.629e-3, init=c(50,0,0), T=500){
#Initiate empty matrix to hold simulation
res = matrix(NA, nrow=500, ncol=3)

#Add column names
dimnames(res) = list(NULL, c("L", "P", "A"))

#store initial conditions in first row
res[1,] = init

for(i in 2:T){
     #Larval equation
     res[i,1] = b*res[i-1,3]*exp(-cel*res[i-1,1]-cea*res[i-1,3])
     #Pupal equation
     res[i,2] = res[i-1,1]*(1-mul)
     #Adult equation
     res[i,3] = res[i-1,2]*exp(-cpa*res[i-1,3])+res[i-1,3]*(1-mua)
     }
return(res)
}

output$plot1 <- renderPlot({
out=lpa(b=input$b, cel=input$cel, cea=input$cea, mul=input$mul,
  cpa=input$cpa, mua=input$mua)
 
 
plot(input$T[1]:input$T[2], out[input$T[1]:input$T[2],3], type = "b", xlab = "Week", 
  ylab = "Abundance", xlim=input$T, ylim=c(0, max(out[input$T[1]:input$T[2],c(1,3)])))
lines(input$T[1]:input$T[2], out[input$T[1]:input$T[2],1], col=2)
 legend("topleft",
        legend=c("A", "L"),
        lty=c(1,1),
        pch=c(1,NA),
         col=c("black", "red"))

   })

output$plot2 <- renderPlot({
out=lpa(b=input$b, cel=input$cel, cea=input$cea, mul=input$mul,
  cpa=input$cpa, mua=input$mua)
 
 
plot(out[input$T[1]:input$T[2],1], out[input$T[1]:input$T[2],3], type = ifelse(input$li==TRUE, "b", "p"), xlab = "L", 
  ylab = "A", xlim=c(0, max(out[input$T[1]:input$T[2],1])), 
  ylim=c(0, max(out[input$T[1]:input$T[2],3])))

   })

output$plot3 <- renderPlot({
out=lpa(b=input$b, cel=input$cel, cea=input$cea, mul=input$mul,
  cpa=input$cpa, mua=input$mua)
 

scatterplot3d(x=out[input$T[1]:input$T[2],1], y=out[input$T[1]:input$T[2],2], 
  z=out[input$T[1]:input$T[2],3], type = ifelse(input$li==TRUE, "b", "p"), xlim=c(0, max(out[input$T[1]:input$T[2],1])),
   ylim=c(0, max(out[input$T[1]:input$T[2],2])),
    zlim=c(0, max(out[input$T[1]:input$T[2],3])),
    xlab="L", ylab="P", zlab="A")

   })


  }
)

#' Launch a shiny-app simulating the Lotka-Volterra competition model
#' @details
#' Launch app for details
#' @examples
#' \dontrun{LVcomp.app}
#' @export
LVcomp.app=shinyApp(ui = fluidPage(
# This creates the User Interface (UI)
  tags$head(tags$style(
    HTML('
         #sidebar1 {
            background-color: #D3D3D3;
        }
    
    #sidebar2 {
      background-color: #ECECEC
    }')
  )),
titlePanel("Lotka-Volterra competition model"),
fluidRow(
column(4, id = "sidebar2",
fluidRow(column(5, id = "sidebar1",
sliderInput("r1", "r1:", 0.3,
               min = 0, max = 1, step=0.01),
sliderInput("a", "alpha12:", 0.3,
               min = 0, max = 1, step=0.01),
sliderInput("K1", "K1:", 150,
               min = 0, max = 300, step=1),
numericInput("N1", "initial N1:", 10,
               min = 0, max = 1000)),
column(5, offset = 1, id = "sidebar1",
sliderInput("r2", "r2:", 0.5,
               min = 0, max = 1, step=0.01),
sliderInput("b", "alpha21:", 0.7,
               min = 0, max = 1, step=0.01),
sliderInput("K2", "K2:", 100,
               min = 0, max = 300, step=1),
numericInput("N2", "initial N2:", 15,
               min = 0, max = 1000)),
column(1)),
fluidRow(
column(6, offset = 3, id = "sidebar1",
numericInput("Tmax", "Tmax:", 100,
               min = 0, max = 200)),
column(3))
),
column(8,  tabsetPanel(
      tabPanel("Simulation", plotOutput("plot1")), 
      tabPanel("Details", 
           withMathJax(
            helpText("MODEL:"),
            helpText("Species 1 $$\\frac{dN_1}{dt} = r_1 N_1 (\\frac{K_1-N_1-\\alpha_{12} N_2}{K_1})$$"),
            helpText("Species 2 $$\\frac{dN_2}{dt} = r_2 N_2 (\\frac{K_2-N_2-\\alpha_{21} N_1}{K_2})$$"),
          helpText("N_1-isocline $$N_2 = \\frac{K_1 - N_1}{\\alpha_{12}}$$"),
          helpText("N_2-isocline $$N_2 = K_2 - \\alpha_{21} N_1$$"),
          helpText("Equilibria $$N_1^* = \\frac{K_1-\\alpha_{12} K_2}{1-\\alpha_{12} \\alpha_{21},
           N_2^* = \\frac{K_2-\\alpha_{21} K_1}{1-\\alpha_{12} \\alpha_{21}}$$")),
           helpText(eval(Attr))
           )
  )
))
),


# This creates the "behind the scenes" code (Server)
server = function(input, output) {
compLV=function(t, y, parameters){
   N1=y[1]
   N2=y[2]

   with(as.list(parameters),{
   dN1 = r1*N1*((K1-N1-a*N2)/K1)
   dN2 = r2*N2*((K2-N2-b*N1)/K2)
   res=c(dN1, dN2)
   list(res)
})
}

output$plot1 <- renderPlot({
N1star=(input$K1-input$a*input$K2)/(1-input$a*input$b)
N2star=(input$K2-input$b*input$K1)/(1-input$a*input$b)


times  = seq(0, input$Tmax, by=0.1)
parms=c(r1=input$r1, r2=input$r2,a=input$a,b=input$b,
K1=input$K1,K2=input$K2)
xstart = c(N1=input$N1, N2=input$N2)

out=ode(y=xstart,
   times=times,
   func=compLV,
   parms=parms)

   out=as.data.frame(out)


par(mfrow=c(1,2))  #This puts two plots side by side each other
plot(times, out$N1, ylab="abundance", xlab="time", type="l",
ylim=range(out[,2:3]))
lines(times, out$N2, col="red")
   legend("topright",
         legend=c("N1", "N2"),
         lty=c(1,1),
          col=c("black", "red"))

plot(NA,xlim=c(0,input$K1*2),ylim=c(0,input$K2*2), xlab="N1", ylab="N2")
fld=flowField(compLV, x.lim=c(0,input$K1*2), y.lim=c(0,input$K2*2), 
parameters=parms, system="two.dim", add=TRUE)

#null clines
curve((input$K1-x)/input$a,col="black",add=TRUE)
curve(input$K2-input$b*x,col="red",add=TRUE)
abline(v=0,col="black")
abline(h=0,col="red")
points(0,0,pch=19)
points(input$K1,0,pch=19)
points(0,input$K2,pch=19)
if(!any(c(N1star, N2star)<0)) points(N1star,N2star,pch=19)
lines(out[,2], out[,3], lwd=2)
    })
   }
)

#' Launch a shiny-app simulating the Lotka-Volterra predation model
#' @details
#' Launch app for details
#' @examples
#' \dontrun{LVpred.app}
#' @export
LVpred.app=shinyApp(
# This creates the User Interface (UI)
ui = fluidPage(
  tags$head(tags$style(
    HTML("
         #sidebar1 {
            background-color: #ECECEC;
        }
    
    #sidebar2 {
      background-color: #ECECEC
    }")
  )),
titlePanel("Lotka-Volterra predation model"),
fluidRow(
column(4, id = "sidebar2",
fluidRow(column(5, id = "sidebar1",
sliderInput("r", "r:", 0.6,
              min = 0, max = 1, step=0.01),
sliderInput("a", "a:", 0.1,
              min = 0, max = 1, step=0.01),
numericInput("N0", "initial N:", 10,
              min = 0, max = 100)),
column(5, offset = 1, id = "sidebar1",
sliderInput("b", "b:", 0.1,
              min = 0, max = 1, step=0.01),
sliderInput("m", "m:", 0.2,
               min = 0, max = 1, step=0.01),
numericInput("P0", "initial P:", 10,
              min = 0, max = 100)),
column(1)),
fluidRow(
column(6, offset = 3, id = "sidebar1",
numericInput("Tmax", "Tmax:", 100,
               min = 0, max = 200)),
column(3))
),
mainPanel(tabsetPanel(
  tabPanel("Simulation", plotOutput("plot1", height = 500)),
     tabPanel("Details", 
     withMathJax(
                helpText("MODEL:"),
          helpText("Prey $$\\frac{dN}{dt} = r N - a N P$$"),
          helpText("Predator $$\\frac{dP}{dt} = b a N P - m P$$"),
          helpText("N-isocline $$P^* = r/a$$"),
          helpText("P-isocline $$N^* = g/b$$"),
          helpText("Equilibria $$N^* = g/b, P^* = r/a$$"))),
           helpText(eval(Attr))
)
)
)
),

# This creates the "behind the scenes" code (Server)
server = function(input, output) {

LV=function(t, y, parameters){
  N=y[1]
  P=y[2]

  with(as.list(parameters),{
  dN = r*N-a*N*P
  dP =b*a*N*P-m*P
  res=c(dN, dP)
  list(res)
})
}

output$plot1 <- renderPlot({
Nstar=input$m/(input$b*input$a)
Pstar=input$r/input$a

times  = seq(0, input$Tmax, by=0.1)
parms=c(r=input$r,a=input$a,b=input$b,m=input$m)
xstart = c(N=input$N0, P=input$P0)
 
out=ode(y=xstart,
  times=times,
  func=LV,
  parms=parms)

  out=as.data.frame(out)

 
par(mfrow=c(1,2))  #This puts two plots side by side each other
plot(times, out$N, ylab="abundance", xlab="time", type="l", ylim=range(out[,1:2]))
lines(times, out$P, col="red")
  legend("right",
        legend=c("N", "P"),
        lty=c(1,1),
         col=c("black", "red"))

plot(out$N, out$P, ylab="predator", xlab="prey", type="l", xlim=range(out[,2]), ylim=range(out[,3]))
abline(h=Pstar, col = "black")
abline(v=Nstar, col = "red")
fld=flowField(LV, x.lim=range(out[,2]), y.lim=range(out[,3]), 
parameters=parms, system="two.dim", add=TRUE)
points(0,0,pch = 1)
points(Nstar,Pstar, pch = 19)
   })
  }
)

#' Launch a shiny-app simulating May's Parasitoid-host Model model
#' @details
#' Launch app for details
#' @examples
#' \dontrun{May.app}
#' @export
May.app=shinyApp(
# This creates the User Interface (UI)
ui = pageWithSidebar(
headerPanel("May's Parasitoid-host Model"),
sidebarPanel(
sliderInput("R", "Growth rate (R):", 1.1,
              min = 1, max = 2, step=.01),
sliderInput("a", "Search efficiency (a):", 0.1,
              min = 0, max = .5),
sliderInput("k", "aggregation (k):", 1.5,
              min = 0.1, max = 3, step=0.1),
numericInput("P0", "Initial parasitoid:", 10,
              min = 1, max = 100),
numericInput("H0", "Initial host:", 20,
              min = 1, max = 100),
numericInput("Tmax", "Tmax:", 100,
              min = 1, max = 500)
),
mainPanel(tabsetPanel(
  tabPanel("Simulation", plotOutput("plot1", height = 500)),
     tabPanel("Details", 
    withMathJax(
         helpText("MODEL:"),
             helpText("Host $$H_t = R H_{t-1} (1 + a P_{t-1})^k$$"),
          helpText("Parasitoid $$P_t = R H_{t-1} (1-(1 + a P_{t-1})^k)$$"),
          helpText("REFERENCE: May RM (1978) Host-parasitoid systems in patchy 
            environments: a phenomenological model. J Anim Ecol 47: 833-843"),
           helpText(eval(Attr))
)
)
)
)
),

# This creates the 'behind the scenes' code (Server)
server = function(input, output) {
 NB = function(R, a, k, T = 100, H0 = 10, P0 = 1){
   #T is length of simulation (number of time-steps)
   #H0 and P0 are initial numbers
   #we provide default parameters except for R and a

   H=rep(NA, T) #host series
   P=rep(NA, T) #parasitoid series

   H[1] = H0 #Initiating the host series
   P[1] = P0 #Initiating the host series

   for(t in 2:T){
     H[t] = R * H[t-1] * (1+ a * P[t-1])^(-k)
     P[t] = R * H[t-1] * (1-(1+ a * P[t-1])^(-k))
     if(P[t-1]==0) break
   } #end of loop

   #the two vectors of results are stored in a "list"
   res= data.frame(H = H, P = P)
 
   #the list is passed out of this function
   return(res)
} #end of function



  output$plot1 <- renderPlot({

    sim= NB(R=input$R, a=input$a, k=input$k, H0=input$H0, P0=input$P0, T=input$Tmax)
    time = 1:input$Tmax

    plot(time, sim$H, type= "b",xlab = "Generations", ylab = "Abundance", 
      ylim = range(sim, na.rm=TRUE))
    points(time, sim$P, type = "b", col = "red")
     legend("topleft",
        legend=c("H", "P"),
        lty=c(1,1),
        pch=c(1,1),
        col=c("black", "red"))
   })
  }
)

#' Launch a shiny-app simulating the Nicholson-Bailey model
#' @details
#' Launch app for details
#' @examples
#' \dontrun{NB.app}
#' @export
NB.app=shinyApp(
# This creates the NB User Interface (UI)
ui = pageWithSidebar(
headerPanel("Nicholson-Bailey Model"),
sidebarPanel(
sliderInput("R", "Growth rate (R):", 1.1,
              min = 1, max = 2, step=.01),
sliderInput("a", "Search efficiency (a):", 0.1,
              min = 0, max = .5),
numericInput("P0", "Initial parasitoid:", 10,
              min = 1, max = 100),
numericInput("H0", "Initial host:", 20,
              min = 1, max = 100),
numericInput("Tmax", "Tmax:", 100,
              min = 1, max = 500)
),
mainPanel(tabsetPanel(
  tabPanel("Simulation", plotOutput("plot1", height = 500)),
  tabPanel("Details",
    withMathJax(
            helpText("MODEL:"),
          helpText("Host $$H_t = R H_{t-1} (1 - \\mbox{exp}(- a P_{t-1}))$$"),
          helpText("Parasitoid $$P_t = R H_{t-1} \\mbox{exp}(- a P_{t-1})$$"),
          helpText("Equilibria $$H^* = \\frac{- \\mbox{log}(R)}{a (R-1)}, 
            P^* = \\frac{- \\mbox{log}(R)}{a}$$"))),
          helpText("REFERENCE: Nicholson AJ, Bailey VA (1935) The balance of animal populations. 
            Proceedings of the Zoological Society of London 3: 551-598"),
           helpText(eval(Attr))
)
)
),


# This creates the 'behind the scenes' code (Server)
server = function(input, output) {
 NB = function(R, a, T = 100, H0 = 10, P0 = 1){
   #T is length of simulation (number of time-steps)
   #H0 and P0 are initial numbers
   #we provide default parameters except for R and a

   H=rep(NA, T) #host series
   P=rep(NA, T) #parasitoid series

   H[1] = H0 #Initiating the host series
   P[1] = P0 #Initiating the host series

   for(t in 2:T){
     H[t] = R * H[t-1] * exp(- a * P[t-1])
     P[t] = R * H[t-1] * (1-exp(- a * P[t-1]))
     if(P[t-1]==0) break
   } #end of loop

   #the two vectors of results are stored in a "list"
   res= data.frame(H = H, P = P)
 
   #the list is passed out of this function
   return(res)
} #end of function



  output$plot1 <- renderPlot({

    sim= NB(R=input$R, a=input$a, H0=input$H0, P0=input$P0, T=input$Tmax)
    time = 1:input$Tmax

    plot(time, sim$H, type= "b",xlab = "Generations", ylab = "Abundance", 
      ylim = range(sim, na.rm=TRUE))
    points(time, sim$P, type = "b", col = "red")
     legend("topleft",
        legend=c("H", "P"),
        lty=c(1,1),
        pch=c(1,1),
        col=c("black", "red"))
   })
  }
)

#' Launch a shiny-app simulating the Ricker model
#' @details
#' Launch app for details
#' @examples
#' \dontrun{Ricker.app}
#' @export
Ricker.app=shinyApp(
# This creates the User Interface (UI)
ui = pageWithSidebar(
headerPanel("Ricker Model"),
sidebarPanel(
sliderInput("r", "Growth rate (r):", 1,
              min = 0, max = 4, step=.1),
sliderInput("K", "Carrying capacity (K):", 100,
              min = 25, max = 200),
numericInput("X0", "Initial number:", 70,
              min = 1, max = 200),
numericInput("Tmax", "Tmax:", 20,
              min = 1, max = 500)
),

mainPanel(tabsetPanel(
  tabPanel("Simulation", plotOutput("plot1", height = 500)),
  tabPanel("Details",
    withMathJax(
                helpText("MODEL:"),
            helpText("$$X_{t+1} = X_t \\mbox{exp}(r (1- X_t/K))$$"),
            helpText("REFERENCE: Ricker WE (1954) Stock and recruitment. 
              Journal of Fishery Research Board Canada 11: 559-623"),
       helpText(eval(Attr))
)
)
)
)
),

# This creates the 'behind the scenes' code (Server)
server = function(input, output) {
 logist = function(r, K, length = 200, X0=70){
  X =  rep(NA, length) #set up the empty vector of the right length
  X[1] = X0 #setting the abundance at time 1 to N0

  for(i in 2:length){ #iteratively updating the growth model.
                    #next abundance is determined by previous abundance
    X[i] = X[i-1]*exp(r*(1-X[i-1]/K))
    }
  return(X) #returning the simulated vector
  }



  output$plot1 <- renderPlot({

    X= logist(r=input$r, K=input$K, length=input$Tmax, X0=input$X0)
    time = 1:input$Tmax
    par(mfrow=c(1,2))
     plot(X, xlab = "time", ylab = "abundance", type="b") # making a time series plot
    curve(x*exp(input$r*(1-x/input$K)),0,input$K*3, xlab = "Xt-1", ylab = "Xt")
abline(a=0, b=1) # adding the 1-to-1 line
points(X[1:(input$Tmax-1)],X[2:input$Tmax], col = "red") # adding the points
# from the simulation to the graph
lines(X[1:(input$Tmax-1)], X[2:input$Tmax], col = "red") # adding the line to connect the points
   })
  }
)


#' Launch a shiny-app simulating the Rosenzweig-MacArthur model
#' @details
#' Launch app for details
#' @examples
#' \dontrun{RM.app}
#' @export
RM.app=shinyApp(
# This creates the User Interface (UI)
ui = fluidPage(
  tags$head(tags$style(
    HTML('
         #sidebar1 {
            background-color: #D3D3D3;
        }
    
    #sidebar2 {
      background-color: #ECECEC
    }')
  )),
titlePanel("Rosenzweig-MacArthur model"),
fluidRow(
column(4, id = "sidebar2",
fluidRow(column(5, id = "sidebar1",
sliderInput("r", "r:", 0.1,
               min = 0, max = 1, step=0.01),
sliderInput("K", "K:", 30,
               min = 0, max = 300, step=1),
sliderInput("a", "a:", 0.2,
               min = 0, max = 1, step=0.01),
numericInput("N", "initial N:", 10,
               min = 0, max = 100)),
column(5, offset = 1, id = "sidebar1",
sliderInput("c", "c:", 20,
               min = 0, max = 100, step=0.1),
sliderInput("b", "b:", 0.1,
               min = 0, max = 1, step=0.01),
sliderInput("g", "g:", 0.05,
               min = 0, max = 1, step=0.01),
numericInput("P", "initial P:", 1,
               min = 0, max = 100)),
column(1)),
fluidRow(
column(6, offset = 3, id = "sidebar1",
numericInput("Tmax", "Tmax:", 1000,
               min = 0, max = 5000)),
column(3))
),
#column(8, plotOutput("plot1", height = 500))
column(8,  tabsetPanel(
      tabPanel("Time", plotOutput("plot1")), 
      tabPanel("Phase plane", plotOutput("plot2")),
      tabPanel("Details", 
           withMathJax(
       helpText("MODEL:"),
            helpText("Prey $$\\frac{dN}{dt} = r N (1-\\frac{N}{K}) - \\frac{a N P}{c+N}$$"),
          helpText("Predator $$\\frac{dP}{dt} = \\frac{b N P}{c+N} - g P$$"),
          helpText("N-isocline $$P^* = (r-rN/K)(c+N)/a$$"),
          helpText("P-isocline $$N^* = gc/(b-g)$$"),
          helpText("Equilibria $$N^* = gc/(b-g), P^* = (r-rN^*/K)(c+N^*)/a$$"),
         helpText("REFERENCE: Rosenzweig ML, MacArthur RH (1963) Graphical representation 
              and stability conditions of predator-prey interactions. Am Nat 97: 209-223"),
           helpText(eval(Attr))
          ))

   
  )
)
)
),

# This creates the "behind the scenes" code (Server)
server = function(input, output){
RM=function(t, y, parameters){
  N=y[1]
  P=y[2]

  r=parameters["r"]
  K=parameters["K"]
  a=parameters["a"]
  c=parameters["c"]
  b=parameters["b"]
  g=parameters["g"]

    dN = r*N*(1-N/K)-a*N*P/(c+N)
    dP = b*N*P/(c+N)-g*P
    res=c(dN,dP)
    list(res)
}

output$plot1 <- renderPlot({

 times  = seq(0, input$Tmax, by=0.1)
 parms=c(r=input$r, K=input$K,a=input$a,
   c=input$c,b=input$b,g=input$g)
 xstart = c(N=input$N, P=input$P)

 out=ode(y=xstart,
    times=times,
    func=RM,
    parms=parms)

    out=as.data.frame(out)

  r=parms["r"]
  K=parms["K"]
  a=parms["a"]
  c=parms["c"]
  b=parms["b"]
  g=parms["g"]


 plot(out$time, out$N, ylab="abundance", xlab="time", type="l", ylim=range(out[,2:3]))
 lines(out$time, out$P, col="red")
    legend("topright",
          legend=c("N", "P"),
          lty=c(1,1),
           col=c("black", "red"))
})

output$plot2 <- renderPlot({

 times  = seq(0, input$Tmax, by=0.1)
 parms=c(r=input$r, K=input$K,a=input$a,
   c=input$c,b=input$b,g=input$g)
 xstart = c(N=input$N, P=input$P)

 out=ode(y=xstart,
    times=times,
    func=RM,
    parms=parms)

    out=as.data.frame(out)

  r=parms["r"]
  K=parms["K"]
  a=parms["a"]
  c=parms["c"]
  b=parms["b"]
  g=parms["g"]

#null clines
plot(out$N, out$P, ylab='predator', xlab='prey', type='l', 
xlim=range(out$N), ylim= range(out$P))
abline(h=0, col = "green")
abline(v=0, col = "red")
curve(r*(1-x/K)*(c+x)/a,from = 0, to = max(c(90, out$N)), col = "green",add=T)
abline(v=g*c/(b-g),col = "red")
fld=flowField(RM, x.lim=range(out$N), y.lim=range(out$P), 
parameters=parms, system="two.dim", add=TRUE)
   legend("topright",
          legend=c("N-iso", "P-iso"),
          lty=c(1,1),
           col=c("green", "red"))

# points(Nstar,Pstar,pch = 1)

})


}
)

#' Launch a shiny-app simulating the seasonal SEIR model
#' @details
#' Launch app for details
#' @examples
#' \dontrun{SEIR.app}
#' @export
SEIR.app=shinyApp(
# This creates the User Interface (UI)
ui = pageWithSidebar(
headerPanel("Seasonally forced SEIR"),
sidebarPanel(
sliderInput("beta0", "Transmission (yr^-1):", 1000,
              min = 0, max = 3000),
sliderInput("beta1", "Seasonality:", 0,
              min = 0, max = 1),
sliderInput("Ip", "Infectious period (days)", 5,
              min = 1, max = 100),
sliderInput("oneoversigma", "Latent period (days):", 8,
              min = 1, max = 100),
sliderInput("mu", "birth rate (per 1000):", 0.02,
              min = 0, max = .1),
sliderInput("T", "Time range:",
                  min = 0, max = 100, value = c(0,20)),
checkboxInput("lg", "un-Log", TRUE)
),
mainPanel(
  tabsetPanel(
      tabPanel("Time", plotOutput("plot1")), 
      tabPanel("Phase plane", plotOutput("plot2")),
       tabPanel("Details", 
           withMathJax(
       helpText("MODEL:"),
            helpText("Susceptible $$\\frac{dS}{dt} = \\mu (N - S) - \\frac{\\beta(t) I S}{N}$$"),
            helpText("Exposed $$\\frac{dE}{dt} = \\frac{\\beta(t) I S}{N} - (\\mu+\\sigma) E$$"),
            helpText("Infectious $$\\frac{dI}{dt} = \\sigma E - (\\mu+\\gamma) I$$"),
           helpText("Removed $$\\frac{dR}{dt} = \\gamma I - \\mu R$$"),
           helpText("Seasonality $$\\beta(t) =  \\beta_0 (1 + \\beta_1 cos(2 \\pi t))$$"),
           helpText("Reproductive ratio $$R_0 =  \\frac{\\sigma}{\\sigma +\\mu} \\frac{1}{\\gamma+\\mu} \\frac{\\beta N}{N}$$"),             
            helpText("REFERENCE: Earn DJD, Rohani P, Bolker BM, Grenfell BT (2000) A simple model for complex dynamical transitions in epidemics.
             Science 287: 667-670"),
           helpText(eval(Attr))
           ))
  
  )
)
),

# This creates the 'behind the scenes' code (Server)
server = function(input, output) {
  seirmod2=function(t, x, params){
  S=x[1]
  E=x[2]
  I=x[3]
  R=x[4]

  mu=params["mu"]
  N=params["N"]
  beta0=params["beta0"]
  beta1=params["beta1"]
  sigma=params["sigma"]
  gamma=params["gamma"]

  dS = mu * (N  - S)  - beta0 * (1+beta1*cos(2*pi*t))* S * I / N
  dE = beta0 * (1+beta1*cos(2*pi * t))* S * I / N - (mu + sigma) * E
  dI = sigma * E - (mu + gamma) * I
  dR = gamma * I - mu * R
  res=c(dS, dE, dI, dR)
  list(res)
} 



  output$plot1 <- renderPlot({

  times  = seq(0, input$T[2], by=1/100)
  paras  = c(mu = input$mu, N = 1, beta0 = input$beta0, beta1 = input$beta1, sigma = 365/input$oneoversigma, gamma = 365/input$Ip)
  xstart = c(S=0.06, E=0, I=0.001, R = 0.939)
  R0 = round(with(as.list(paras), sigma/(sigma+mu)*beta0/(gamma+mu)), 1)
 
out=ode(y=xstart,
  times=times,
  func=seirmod2,
  parms=paras)

  out=as.data.frame(out)

  sel=out$time>input$T[1]&out$time<input$T[2]

par(mar = c(5,5,2,5))
#lg=ifelse(input$lg==TRUE, "y", "")
plot(x=out$time[sel], y=out$I[sel], ylab="fraction", xlab="time", type="l",
ylim=range(out[sel,-c(1,2, 5)]), xlim=c(input$T[1], input$T[2]), log=ifelse(input$lg==TRUE, "y", ""), col="red")
 lines(x=out$time, y=out$E, col="blue")
title(paste("R0=", R0))
# lines(x=out$time, y=out$S, col="green")
par(new=T)
plot(x=out$time, y=out$S, type="l", col="green", axes=FALSE, xlab=NA, ylab=NA, 
    ylim=range(out[sel,2]), xlim=c(input$T[1], input$T[2]), log=ifelse(input$lg==TRUE, "y", ""))
axis(side = 4, col="green")
mtext(side = 4, line = 4, "S", col="green")
  legend("right",
        legend=c("I", "E", "S"),
        lty=c(1,1,1),
         col=c("red", "blue", "green"))
   })
  
output$plot2 <- renderPlot({
  times  = seq(0, input$T[2], by=1/100)
  paras  = c(mu = input$mu, N = 1, beta0 = input$beta0, beta1 = input$beta1, sigma = 365/input$oneoversigma, gamma = 365/input$Ip)
  xstart = c(S=0.06, E=0, I=0.001, R = 0.939)
  R0 = with(as.list(paras), sigma/(sigma+mu)*beta0/(gamma+mu))
 
  out=ode(y=xstart,
  times=times,
  func=seirmod2,
  parms=paras)

  out=as.data.frame(out)

  sel=out$time>input$T[1]&out$time<input$T[2]

  plot(out$S[sel], out$I[sel], log=ifelse(input$lg==TRUE, "xy", ""), type="l", xlab="fraction susceptible", ylab="fraction infected")
  abline(v=1/R0, col="green")
  curve(paras["mu"]*(1-x)/(paras["beta0"]*x), min(out$S), max(out$S), add=TRUE, col="red")
    legend("topright",
        legend=c("I-socline", "S-isocline"),
        lty=c(1,1),
         col=c("red", "green"))
 
  })

  }
)

#' Launch a shiny-app simulating the SEIRS model
#' @details
#' Launch app for details
#' @examples
#' \dontrun{SEIRS.app}
#' @export
SEIRS.app=shinyApp(
# This creates the User Interface (UI)
ui <- pageWithSidebar(
headerPanel("SEIRS periodicity"),
sidebarPanel(
sliderInput("beta", "Transmission (yr^-1):", 500,
              min = 0, max = 3000),
sliderInput("oneoveromega", "Immune duration (years):", 4,
              min = 0, max = 100),
sliderInput("Ip", "Infectious period (days)", 5,
              min = 1, max = 100),
sliderInput("oneoversigma", "Latent period (days):", 8,
              min = 1, max = 100),
sliderInput("oneovermu", "Life expectancy (years):", 10,
              min = 1, max = 100),
sliderInput("T", "Time range:",
                  min = 0, max = 100, value = c(0,20)),
checkboxInput("lg", "un-Log", TRUE)
),

mainPanel(
  tabsetPanel(
      tabPanel("Time", plotOutput("plot1")), 
      tabPanel("Phase plane", plotOutput("plot2")),
       tabPanel("Equations", 
           withMathJax(
            helpText("Susceptible $$\\frac{dS}{dt} = \\mu (N - S) - \\frac{\\beta I S}{N} + \\omega R$$"),
            helpText("Exposed $$\\frac{dE}{dt} = \\frac{\\beta I S}{N} - (\\mu+\\sigma) E$$"),
            helpText("Infecitous $$\\frac{dI}{dt} = \\sigma E - (\\mu+\\gamma) I$$"),
           helpText("Removed $$\\frac{dR}{dt} = \\gamma I - \\mu R - \\omega R$$"),
           helpText("Reproductive ratio $$R_0 =  \\frac{\\sigma}{\\sigma +\\mu} \\frac{1}{\\gamma+\\mu} \\frac{\\beta N}{N}$$"),
           helpText(eval(Attr))             
           ))
  
  )
)
),

# This creates the 'behind the scenes' code (Server)
server <- function(input, output) {
seirsmod=function(t, x, params){
  S=x[1]
  E=x[2]
  I=x[3]
  R=x[4]

  mu=params["mu"]
  beta=params["beta"]
  omega=params["omega"]
  sigma=params["sigma"]
  gamma=params["gamma"]

  dS = mu * (1  - S)  - beta * S * I + omega * R
  dE = beta * S * I - (mu + sigma) * E
  dI = sigma * E - (mu + gamma) * I
  dR = gamma * I - mu * R - omega * R
  res=c(dS, dE, dI, dR)
  list(res)
} 

  output$plot1 <- renderPlot({

  times  = seq(0, input$T[2], by=1/100)
  paras  = c(mu = 1/input$oneovermu, beta =  input$beta, sigma = 365/input$oneoversigma, gamma = 365/input$Ip, omega=1/input$oneoveromega)
  xstart = c(S=0.539, E=0, I=0.001, R = 0.46)

  R0=with(as.list(paras),{
    beta*sigma/((mu+sigma)*(mu+gamma))
  })


Sstar=1/R0
Istar=paras["mu"]*(1-Sstar)/(paras["beta"]*Sstar - (paras["omega"]*paras["gamma"])/(paras["mu"]+paras["omega"]))
Estar=(paras["mu"]+paras["gamma"])*Istar/paras["sigma"]
Rstar=paras["gamma"]*Istar/(paras["mu"]+paras["omega"])

star=as.list(c(S=Sstar, E=Estar, I=Istar, R=Rstar, paras))
names(star)[1:4]=c("S", "E", "I", "R")

fns=list(quote(mu * (1  - S)  - beta * S * I  + omega * R), quote(beta * S * I - (mu + sigma) * E), quote(sigma * E - (mu + gamma) * I), quote(gamma * I - mu * R - omega * R))

aa1=as.vector(sapply(fns, D, "S"))
aa2=as.vector(sapply(fns, D, "E"))
aa3=as.vector(sapply(fns, D, "I"))
aa4=as.vector(sapply(fns, D, "R"))

JJ=matrix(c(sapply(aa1, eval, star), sapply(aa2, eval, star),sapply(aa3, eval, star),sapply(aa4, eval, star)), ncol=4)

EE=eigen(JJ)$values
WW=which.max(Im(EE))
rp=2*pi/Im(EE[WW])


out=ode(y=xstart,
  times=times,
  func=seirsmod,
  parms=paras)

  out=as.data.frame(out)

  sel=out$time>input$T[1]&out$time<input$T[2]

par(mar = c(5,5,2,5))
plot(x=out$time[sel], y=out$I[sel], ylab="fraction", xlab="time", type="l",
ylim=range(out[sel,-c(1,2, 5)]), xlim=c(input$T[1], input$T[2]), log=ifelse(input$lg==TRUE, "y", ""), col="red")
 lines(x=out$time, y=out$E, col="blue")
 title(paste("R0=", round(R0, 1), ", Period=", round(rp,2)))

par(new=T)
plot(x=out$time, y=out$S, type="l", col="green", axes=FALSE, xlab=NA, ylab=NA, 
    ylim=range(out[sel,2]), xlim=c(input$T[1], input$T[2]), log=ifelse(input$lg==TRUE, "y", ""))
axis(side = 4, col="green")
mtext(side = 4, line = 4, "S", col="green")
  legend("right",
        legend=c("I", "E", "S"),
        lty=c(1,1,1),
         col=c("red", "blue", "green"))
   })
  
output$plot2 <- renderPlot({
  times  = seq(0, input$T[2], by=1/100)
  paras  = c(mu = 1/input$oneovermu, beta =  input$beta, sigma = 365/input$oneoversigma, gamma = 365/input$Ip, omega=1/input$oneoveromega)
  xstart = c(S=0.539, E=0, I=0.001, R = 0.46)

  R0=with(as.list(paras),{
    beta*sigma/((mu+sigma)*(mu+gamma))
  })
 

Sstar=1/R0
Istar=paras["mu"]*(1-Sstar)/(paras["beta"]*Sstar - (paras["omega"]*paras["gamma"])/(paras["mu"]+paras["omega"]))
Estar=(paras["mu"]+paras["gamma"])*Istar/paras["sigma"]
Rstar=paras["gamma"]*Istar/(paras["mu"]+paras["omega"])

star=as.list(c(S=Sstar, E=Estar, I=Istar, R=Rstar, paras))
names(star)[1:4]=c("S", "E", "I", "R")

fns=list(quote(mu * (1  - S)  - beta * S * I  + omega * R), quote(beta * S * I - (mu + sigma) * E), quote(sigma * E - (mu + gamma) * I), quote(gamma * I - mu * R - omega * R))

aa1=as.vector(sapply(fns, D, "S"))
aa2=as.vector(sapply(fns, D, "E"))
aa3=as.vector(sapply(fns, D, "I"))
aa4=as.vector(sapply(fns, D, "R"))

JJ=matrix(c(sapply(aa1, eval, star), sapply(aa2, eval, star),sapply(aa3, eval, star),sapply(aa4, eval, star)), ncol=4)

EE=eigen(JJ)$values
WW=which.max(Im(EE))
rp=2*pi/Im(EE[WW])


  out=ode(y=xstart,
  times=times,
  func=seirsmod,
  parms=paras)

  out=as.data.frame(out)

  sel=out$time>input$T[1]&out$time<input$T[2]

  plot(out$S[sel], out$I[sel], log=ifelse(input$lg==TRUE, "xy", ""), type="l", xlab="fraction susceptible", ylab="fraction infected")
  title(paste("R0=", round(R0, 1), ", Period=", round(rp,2)))
  abline(v=1/R0, col="green")
  curve(paras["mu"]*(1-x)/(paras["beta"]*x - (paras["omega"]*paras["gamma"])/(paras["mu"]+paras["omega"])), min(out$S), max(out$S), add=TRUE, col="red")
    legend("topright",
        legend=c("I-socline", "S-isocline"),
        lty=c(1,1),
         col=c("red", "green"))
 
  })

  }
)

#' Launch a shiny-app simulating the SIR model
#' @details
#' Launch app for details
#' @examples
#' \dontrun{SIR.app}
#' @export
SIR.app=shinyApp(
# This creates the User Interface (UI)
ui <- pageWithSidebar(
headerPanel("The SIR model"),
sidebarPanel(
sliderInput("beta", "Transmission (yr^-1):", 300,
              min = 0, max = 1000),
sliderInput("infper", "Infectious period (days)", 5,
              min = 1, max = 100),
sliderInput("mu", "birth rate (/year):", 5,
              min = 0, max = 100),
sliderInput("T", "Time range:",
                  min = 0, max = 1, value = c(0,1))
),
mainPanel(
  tabsetPanel(
      tabPanel("Time", plotOutput("plot1")), 
      tabPanel("Phase plane", plotOutput("plot2", height = 500)),
      tabPanel("Equations", 
           withMathJax(
            helpText("Susceptible $$\\frac{dS}{dt} = \\mu (N - S) - \\frac{\\beta I S}{N}$$"),
            helpText("Infecitous $$\\frac{dI}{dt} = \\frac{\\beta I S}{N} - (\\mu+\\sigma) I$$"),
           helpText("Removed $$\\frac{dR}{dt} = \\gamma I - \\mu R$$"),
           helpText("Reproductive ratio $$R_0 =  \\frac{1}{\\gamma+\\mu} \\frac{\\beta N}{N}$$")             
           ))
  
)   
  )
 )
,


# This creates the 'behind the scenes' code (Server)
server <- function(input, output) {
  sirmod=function(t, x, parms){
    S=x[1]
    I=x[2]
    R=x[3]

    beta=parms["beta"]
    mu=parms["mu"]
    gamma=parms["gamma"]
    N=parms["N"]

    dS = mu * (N  - S)  - beta * S * I / N
    dI = beta * S * I / N - (mu + gamma) * I
    dR = gamma * I - mu * R
    res=c(dS, dI, dR)
    list(res)
  }

  output$plot1 <- renderPlot({
  times  = seq(0, input$T[2], by=1/1000)
  parms  = c(mu = input$mu, N = 1, beta =  input$beta, gamma =
    365/input$infper)
  start = c(S=0.999, I=0.001, R = 0)
  R0 = round(with(as.list(parms), beta/(gamma+mu)), 1)

  AA=with(as.list(parms), 1/(mu*(R0-1)))
  GG=with(as.list(parms), 1/(mu+gamma))
  rp=round(2*pi*sqrt(AA*GG),2)

  out=ode(y=start,
  times=times,
  func=sirmod,
  parms=parms)

  out=as.data.frame(out)

  sel=out$time>input$T[1]&out$time<input$T[2]

  plot(x=out$time[sel], y=out$S[sel], ylab="fraction", xlab="time", type="l",
  ylim=range(out[sel,-c(1,4)]))
  title(paste("R0=", R0, "Period=", rp))
  lines(x=out$time[sel], y=out$I[sel], col="red")
  lines(x=out$time[sel], y=out$R[sel], col="green")
  legend("right",
        legend=c("S", "I", "R"),
        lty=c(1,1,1),
         col=c("black", "red", "green"))
   })

  output$plot2 <- renderPlot({
  times  = seq(0, input$T[2], by=1/1000)
  parms  = c(mu = input$mu, N = 1, beta =  input$beta, gamma =
    365/input$infper)
  start = c(S=0.999, I=0.001, R = 0)
  R0 = with(as.list(parms), beta/(gamma+mu))

AA=with(as.list(parms), 1/(mu*(R0-1)))
  GG=with(as.list(parms), 1/(mu+gamma))
  rp=round(2*pi*sqrt(AA*GG),2)
 
simod=function(t, y, parameters){
   S=y[1]
   I=y[2]

   beta=parameters["beta"]
   mu=parameters["mu"]
   gamma=parameters["gamma"]
   N=parameters["N"]
   
   dS = mu * (N  - S)  - beta * S * I / N
   dI = beta * S * I / N - (mu + gamma) * I
   res=c(dS, dI)
   list(res)
 }

  out=ode(y=start[-3],
  times=times,
  func=simod,
  parms=parms)

  out=as.data.frame(out)

  plot(x=out$S, y=out$I, xlab="Fraction suceptible", ylab="Fraction infected", type="l")
  title(paste("R0=", round(R0, 1), "Period=", rp))
  
  require(phaseR)
fld=flowField(simod, x.lim=range(out$S), y.lim=range(out$I), 
parameters=parms, system="two.dim", add=TRUE,
ylab="I", xlab="S")
#cli=nullclines(simod, x.lim=c(0,.4), y.lim=c(0,.01), 
#parameters=parms, system="two.dim", add=TRUE, points=201)


  abline(v=1/R0, col="green")
  curve(parms["mu"]*(1-x)/(parms["beta"]*x), min(out$S), max(out$S), add=TRUE, col="red")
    legend("topright",
        legend=c("I-socline", "S-isocline"),
        lty=c(1,1),
         col=c("red", "green"))

   })
  }
)

#' Launch a shiny-app simulating TSIR model
#' @details
#' Launch app for details
#' @examples
#' \dontrun{TSIR.app}
#' @export
TSIR.app=shinyApp(
# This creates the User Interface (UI)
ui = pageWithSidebar(
headerPanel("Simulating with TSIR"),
sidebarPanel(
sliderInput("alpha", "alpha:", 0.97,
              min = 0.8, max = 1),
sliderInput("beta", "beta:", 25,
              min = 0, max = 100),
sliderInput("B", "Births (B):", 2300,
              min = 0, max = 5000),
sliderInput("S0", "Fraction S:", 0.06,
              min = 0, max = 1),
numericInput("sdbeta", "sdbeta:", 3,
              min = 0, max = 10),
numericInput("N", "Popsize:", 3000000,
              min = 0, max = 100),
numericInput("I0", "initial I:", 100,
              min = 1, max = 100),
numericInput("IT", "Iterations:", 520,
              min = 1, max = 1000)
),
mainPanel(
  tabsetPanel(
      tabPanel("Simulation", plotOutput("plot1")), 
      tabPanel("Transfer function", plotOutput("plot2")),
      tabPanel("Details", 
           withMathJax(
            helpText("MODEL:"),
            helpText("Susceptible $$S_{t+1} = S_t - I_{t+1} + B$$"),
            helpText("Expected Inf $$\\lambda_{t+1} = \\frac{\\beta_t I_t^\\alpha S_t}{N}$$"),
           helpText("Infected $$I_{t+1} \\sim \\mbox{Poisson}(\\lambda_{t+1})$$"),
           helpText("Transmission $$\\beta_t \\sim \\mbox{Norm}(\\beta, \\mbox{sd}\\beta^2)$$"),             
            helpText("The transfer function is $$T(\\omega) = (\\vec{I}-e^{-\\imath \\omega} \\vec{J})^{-1} \\cdot \\vec{A}$$")
          )))   
  )),

# This creates the 'behind the scenes' code (Server)
server = function(input, output) {
  SimTsir=function(alpha, B, beta, sdbeta,
    S0, I0, IT, N){
    lambda = rep(NA, IT)
    I = rep(NA, IT)
    S = rep(NA, IT)
    I[1] = I0
    lambda[1] = I0
    S[1] = S0*N
    for(i in 2:IT) {
        lambda[i] = rnorm(1, mean=beta, sd=sdbeta) * I[i - 1]^alpha * S[i - 1] /N
        if(lambda[i]<0) {lambda[i]=0}
        I[i] = rpois(1, lambda[i])
        S[i] = S[i - 1] + B - I[i]
    }
    list(I = I, S = S)
}

  output$plot1 <- renderPlot({
  out = SimTsir(alpha=input$alpha, B=input$B, beta=input$beta, sdbeta=input$sdbeta, 
  S0=input$S0, I0=input$I0, IT=input$IT, N=input$N)
par(mfrow=c(1,2))  #This puts two plots side by side each other
plot(out$I, ylab="infected", xlab="time", type="b")
plot(out$S, out$I, ylab="infected", xlab="susceptible", type="b")
})

  output$plot2 <- renderPlot({
  out = SimTsir(alpha=input$alpha, B=input$B, beta=input$beta, sdbeta=input$sdbeta, 
  S0=input$S0, I0=input$I0, IT=input$IT, N=input$N)
Seq=expression(S-beta*S*I^alpha/N+B)
Ieq=expression(beta*S*I^alpha/N)
j11=D(Seq, "S")
j12=D(Seq, "I")
j21=D(Ieq, "S")
j22=D(Ieq, "I")
jj=c(j11, j12, j21,j22)

a1=D(Seq, "beta")
a2=D(Ieq, "beta")
aa=c(a1, a2)

paras  = c(B = input$B, beta =  input$beta, alpha = input$alpha, N=input$N)
eqs=sapply(c(quote(B^(1-alpha)*N/beta), quote(B)), eval, as.list(paras))

J=matrix(sapply(jj, eval, as.list(c(paras, c(S=eqs[1], I=eqs[2])))), 2, byrow=TRUE)
evs=eigen(J)$values
rp=2*pi/atan2(Im(evs[1]), Re(evs[1]))

A=matrix(sapply(aa, eval, as.list(c(paras, c(S=eqs[1], I=eqs[2])))), 2, byrow=TRUE)
Id=matrix(c(1,0,0,1),ncol=2)
wseq=seq(0,pi,length=500)

Fr=vector("list",500)  #set up empty list of matrices
# loop to fill those matrices with fourier trasform for the 500 values of w
for(i in 1:500){ 
Fr[[i]]=matrix(solve(Id-exp(1i*wseq[i])*J)%*%A,ncol=1)  #solve gives inverse
}

PS=matrix(NA,ncol=2,nrow=500,dimnames=list(1:500, c("S","I"))) 
#power spectra from real and imaginary parts of Fourier transform
for(i in 1:500){
PS[i,]=sqrt(Re(Fr[[i]])^2+Im(Fr[[i]])^2)
}

sfit=spectrum(out$I[-c(1:104)])

require(polspline)
sfit2=lspec(out$I[-c(1:104)])
plot(wseq, PS[,2], type="l", 
	xlab="frequency (in radians)", ylab="amplitude", xlim=c(0,0.6))
title("Simulated spectrum (periodogram and log-spline) \n and T-fn prediction")
lines(pi*sfit$freq/0.5, max(PS[,2])*sfit$spec/max(sfit$spec), col=2)
par(new=TRUE)
plot(sfit2, col=3, xlim=c(0,0.6), axes=FALSE)
legend("topright", c("T-function", "periodogram", "log-spline"), lty=c(1,1,1), col=c(1,2,3))
   })
  }
)

