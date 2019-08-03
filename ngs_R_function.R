##################################################################################################################################
# Data 2019-08-03
# Author Howard MENG
# E-mail meng_howard@126.com
# useful R function  
##################################################################################################################################

#------------------------------------------------------------------------------------------>>>>>>>
# contents
#------------------------------------------------------------------------------------------>>>>>>>
# plot.TAD function for plot Hi-C TAD
# plot.matrix easy to plot heatmap 

#------------------------------------------------------------------------------------------>>>>>>>
# function part
#------------------------------------------------------------------------------------------>>>>>>>
plot.TAD <- function(input.matrix,col.min = "red",col.max = "red",col.boundary = 0,ylim=NULL,minBound=0,maxBound=0.95,mat.upper=T,mat.part=T){
  
  # 参数说明
  ## input.matrix: 输入矩阵
  ## col.min="red": 热图颜色，支持字符或者是RGB参数za
  ## col.max="red": 热图颜色，支持字符或者是RGB参数
  ## col.boundary=0: 两种颜色的分界线，也即白色对应的值
  ## ylim=NULL: 热图的y轴范围，1个bin的长度是1，方便截短
  ## minBound=0: 使用分位数来修正最小值
  ## maxBound=0.95: 使用分位数修正最大值，默认matrix中的95%分位数为最大值
  ## mat.upper=T: 如果为TRUE画上半部分，为FALSE画下半部分; 只有在mat.part=T时有效
  ## mat.part=T: 如果为TRUE则画三角，如果为FLASE则画倒置的matrix
  input.matrix = as.matrix(input.matrix)
  mat.nrow = dim(input.matrix)[1]
  
  ## color type 1：single color
  ## color type 2：blue white red black
  
  ## 生成矩阵需要的颜色
  mat.lower.tri = as.vector(input.matrix[lower.tri(input.matrix,diag = F)])
  mat.quantile = quantile(as.vector(input.matrix[lower.tri(input.matrix,diag = T)]),prob=c(minBound,maxBound))
  
  if(col.boundary>mat.quantile[2] | col.boundary<mat.quantile[1]){
    print("Error! col.boundary have to smaller than matrix maxBound value!")
    return(NULL)
  }
  
  ## fix value,too large or too small
  mat.lower.tri[mat.lower.tri < mat.quantile[1]] = mat.quantile[1]
  mat.lower.tri[mat.lower.tri > mat.quantile[2]] = mat.quantile[2]
  mat.diagnal.value = diag(input.matrix)
  mat.diagnal.value[mat.diagnal.value < mat.quantile[1]] = mat.quantile[1]
  mat.diagnal.value[mat.diagnal.value > mat.quantile[2]] = mat.quantile[2]
  
  ## create color vector
  mat.lower.tri.col_alpha = rep(0,length(mat.lower.tri))
  mat.diagnal.col_alpha = rep(0,length(mat.diagnal.value))
  
  mat.lower.tri.col_alpha[mat.lower.tri>=col.boundary] = ceiling((mat.lower.tri[mat.lower.tri>=col.boundary] - col.boundary) / (mat.quantile[2] - col.boundary) * 255)
  mat.lower.tri.col_alpha[mat.lower.tri<col.boundary] = ceiling((col.boundary - mat.lower.tri[mat.lower.tri<col.boundary]) / (col.boundary - mat.quantile[1]) * 255)
  mat.diagnal.col_alpha[mat.diagnal.value>=col.boundary] = ceiling((mat.diagnal.value[mat.diagnal.value>=col.boundary] - col.boundary) / (mat.quantile[2] - col.boundary) * 255)
  mat.diagnal.col_alpha[mat.diagnal.value<col.boundary] = ceiling((col.boundary - mat.diagnal.value[mat.diagnal.value<col.boundary]) / (col.boundary - mat.quantile[1]) * 255)
  
  mat.lower.tri.color = rep("##FFFFFF",length(mat.lower.tri))
  mat.diagnal.color = rep("##FFFFFF",length(mat.diagnal.value))
  
  mat.lower.tri.color[mat.lower.tri>=col.boundary] = rgb(t(col2rgb(col.max)),alpha = mat.lower.tri.col_alpha[mat.lower.tri>=col.boundary],maxColorValue = 255)
  mat.lower.tri.color[mat.lower.tri<col.boundary] = rgb(t(col2rgb(col.min)),alpha = mat.lower.tri.col_alpha[mat.lower.tri<col.boundary],maxColorValue = 255)
  mat.diagnal.color[mat.diagnal.value>=col.boundary] = rgb(t(col2rgb(col.max)),alpha = mat.diagnal.col_alpha[mat.diagnal.value>=col.boundary],maxColorValue = 255)
  mat.diagnal.color[mat.diagnal.value<col.boundary] = rgb(t(col2rgb(col.min)),alpha = mat.diagnal.col_alpha[mat.diagnal.value<col.boundary],maxColorValue = 255)
  
  # 根据矩阵的行列，生成x1,y1的坐标
  x1.part = c(1:(mat.nrow -1))
  y1.part = c(1:(mat.nrow -1))
  x1 = x1.part
  y1 = y1.part
  for(i in c(1:(mat.nrow-1))){
    x1.part = x1.part[-length(x1.part)] + 2
    y1.part = y1.part[-length(y1.part)]
    
    x1 = c(x1,x1.part)
    y1 = c(y1,y1.part)
  }
  
  # 对角线的x1
  x1.diagonal = c(0:(mat.nrow-1)) * 2
  y1.diagonal = rep(0,mat.nrow)
  
  if(mat.part==T){
    # plot matrix upper OR lower only.
    x1 = c(x1,x1.diagonal)
    y1 = c(y1,y1.diagonal)  
    # 颜色向量也不同
    mat.color = c(mat.lower.tri.color,mat.diagnal.color)
  }else if(mat.part==F){
    x1 = c(x1,x1.diagonal,x1)
    y1 = c(y1,y1.diagonal,-y1)  
    mat.color = c(mat.lower.tri.color,mat.diagnal.color,mat.lower.tri.color)
  }
  
  
  # 根据x1 y1 生成剩余坐标
  x2 = x1 + 1
  x3 = x2 + 1
  x4 = x2
  y2 = y1 + 1
  y3 = y1
  y4 = y1 - 1
  # 需要按照x1,x2,x3,x4,NA的格式进行生成x vector否则会把所有的点连在一起
  NA_vector = rep(NA,length(x1))
  # 生成按照x11,x12,x13,x14,NA,x21,x22,x23,x24,NA...排列的x与y
  x_matrix = matrix(c(x1,x2,x3,x4,NA_vector),ncol = 5)
  y_matrix = matrix(c(y1,y2,y3,y4,NA_vector),ncol = 5)
  x = as.vector(t(x_matrix))
  y = as.vector(t(y_matrix))
  
  # 绘图部分
  ## 如果mat.upper 为F 则画倒置的图像
  ## 设置画布大小
  x.point = c(0,mat.nrow*2)
  if(mat.part){
    if(! mat.upper){ y = -y ; y.point = c(-mat.nrow,1)}else{ y.point = c(-1,mat.nrow) }
  }else{ 
    y.point = c(-mat.nrow,mat.nrow)
  }
  ## 设置ylim
  if(is.null(ylim)){ ylim = y.point }
  
  ## 生成画布
  plot(x.point,y.point,ylim=ylim,type="n",frame.plot =F,xaxt="n",yaxt="n",xlab="",ylab="")
  ## 画三角矩阵
  polygon(x,y,col = mat.color,border = F)
}


plot.matrix <- function(mat,bound.min=0,bound.max=1,mat.lim=NULL,color_type=1,col.min = "red",col.max = "red",col.boundary = NULL){
  # mat = hic matrix
  # min_bound min quantile to miss data
  # max_bound max quantile to miss data
  
  mat <- as.matrix(mat)
  
  #matrix info calculate
  row_num <- dim(mat)[1]
  col_num <- dim(mat)[2]
  
  #matrix rects' coordinate
  x1 <- rep(c(0:(col_num-1)),each=row_num)
  x2 <- x1 + 1
  y1 <- rep(c(0:(-row_num+1)),col_num)
  y2 <- y1 -1
  
  #matrix colour vector
  if(color_type == 1){
    ## 生成矩阵需要的颜色
    ###确定 matrix的上下界
    mat.quantile = quantile(as.vector(mat),prob=c(bound.min,bound.max))
    
    if(is.null(mat.lim)){
      mat.lim = mat.quantile
    }else{
      mat.lim = c(min(c(mat.quantile,mat.lim)),max(c(mat.quantile,mat.lim)))
    }
    
    if(is.null(col.boundary)){
      col.boundary = mat.lim[1]
    }
    
    if(col.boundary>mat.lim[2]){
      print("Error! col.boundary have to smaller than matrix mat.lim[2] value!")
      return(NULL)
    }else if(col.boundary<mat.lim[1]){
      print("Error! col.boundary have to larger than matrix mat.lim[1] value!")
      print(mat.lim)
      print(col.boundary)
      return(NULL)
    }
    
    ## fix value,too large or too small
    mat[mat < mat.lim[1]] = mat.lim[1]
    mat[mat > mat.lim[2]] = mat.lim[2]
    
    ## create color vector
    mat.col_alpha = rep(0,length(as.vector(mat)))
    mat.col_alpha[mat >= col.boundary] = ceiling((mat[mat >= col.boundary] - col.boundary) / (mat.lim[2] - col.boundary) * 255)
    mat.col_alpha[mat < col.boundary] = ceiling((col.boundary - mat[mat < col.boundary]) / (col.boundary - mat.lim[1]) * 255)
    
    mat.color = rep("#FFFFFF",length(as.vector(mat)))
    mat.color[mat>=col.boundary] = rgb(t(col2rgb(col.max)),alpha = mat.col_alpha[mat>=col.boundary],maxColorValue = 255)
    mat.color[mat< col.boundary] = rgb(t(col2rgb(col.min)),alpha = mat.col_alpha[mat< col.boundary],maxColorValue = 255)
    
  }else if(color_type == 2){
    # input a col_list 
  }
  
  # plot the final matrix
  plot(x=c(0,col_num),y=c(-row_num,0),type="n",frame.plot = F,xaxt="n",yaxt="n",cex.main = 2,xlab="",ylab="")
  rect(x1,y1,x2,y2,col = mat.color,border = NA)
}



#------------------------------------------------------------------------------------------>>>>>>>
# MIT License
#------------------------------------------------------------------------------------------>>>>>>>
# Copyright (c) 2019 MENG Howard
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


# 2019-08-03 By MENG

