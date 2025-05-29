## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(bases)
data(mtcars)

m = lm(mpg ~ b_rff(cyl, disp, hp, wt, p = 10), mtcars)

## -----------------------------------------------------------------------------
all.equal(predict(m), predict(m, newdata = mtcars))
all.equal(predict(m, newdata = mtcars[5:10, ]), 
          predict(m, newdata = mtcars[5:10, ]))

## -----------------------------------------------------------------------------
B = with(mtcars, b_rff(cyl, disp, hp, wt, p = 10))

all.equal(B, predict(B))
all.equal(B, predict(B, newdata = mtcars), check.attributes = FALSE)
nrow(predict(B, newdata = mtcars[1:3, ]))

