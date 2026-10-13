set.seed(1234)

## simple and bipartite samplers (they inherit the plot methods of SimpleSBM / BipartiteSBM)
simpleSampler <- sampleSimpleSBM(nbNodes = 40, blockProp = c(.4, .6),
                                 connectParam = list(mean = matrix(c(.8, .1, .1, .4), 2, 2)))
bipartiteSampler <- sampleBipartiteSBM(nbNodes = c(30, 50),
                                       blockProp = list(row = c(.5, .5), col = c(.2, .3, .5)),
                                       connectParam = list(mean = matrix(c(.9, .4, .1, .2, .6, .05), 2, 3, byrow = TRUE)))

test_that("matrix plots work for all types", {
  for (mySampler in list(simpleSampler, bipartiteSampler)) {
    for (type in c("data", "expected", "graphon")) {
      expect_s3_class(plot(mySampler, type = type), "ggplot")
    }
    expect_s3_class(plot(mySampler, type = "data", ordered = FALSE), "ggplot")
    expect_s3_class(plot(mySampler, type = "data", plotOptions = list(blockBackground = FALSE)), "ggplot")
  }
})

test_that("the block background is drawn only with a clustering", {
  ## with a clustering, the background is a first layer of rectangles (GeomTile inherits from GeomRect)
  geoms <- function(g) unname(sapply(g$layers, function(l) class(l$geom)[1]))
  expect_equal(geoms(plot(simpleSampler, type = "data"))[1], "GeomRect")
  expect_false("GeomRect" %in% geoms(plot(simpleSampler, type = "data", ordered = FALSE)))
  expect_false("GeomRect" %in% geoms(plot(simpleSampler, type = "expected")))
})

test_that("plotMyGraphon works for simple and bipartite SBM", {
  g <- plotMyGraphon(simpleSampler$connectParam, simpleSampler$blockProp, simpleSampler$nbNodes, dimLabels = "node")
  expect_s3_class(g, "ggplot")
  ## one rectangle per pair of blocks
  expect_equal(nrow(g$data), 4)
  expect_equal(sort(g$data$value), sort(c(simpleSampler$connectParam$mean)))

  g <- plotMyGraphon(bipartiteSampler$connectParam$mean, bipartiteSampler$blockProp, bipartiteSampler$nbNodes,
                     dimLabels = c(row = "row", col = "col"), plotOptions = list(legend = TRUE, limits = c(0, 1)))
  expect_s3_class(g, "ggplot")
  expect_equal(nrow(g$data), 6)
  ## the blocks cover the whole network
  expect_equal(max(g$data$xmax), unname(bipartiteSampler$nbNodes[2]))
  expect_equal(min(g$data$ymin), -unname(bipartiteSampler$nbNodes[1]))
})

test_that("plotMyGraphon checks the dimensions", {
  expect_error(plotMyGraphon(matrix(.5, 2, 2), blockProp = c(.2, .3, .5)),
               "do not match")
  expect_error(plotMyGraphon(matrix(.5, 2, 3), blockProp = list(c(.5, .5), c(.5, .5))),
               "do not match")
})
