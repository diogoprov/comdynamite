web <- food_web_dynamics(seq(0, 2000, by = 20), init = c(R = 50, N = 0.1, P = 0),
                         traits = c(R = 10, N = 4, P = 6), dt = 0.5)
tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")

test_that("food_web_rates() gives each tip the speed of its level", {
  r <- food_web_rates(web, tree, c(a = "R", b = "N", c = "R"), sigma2 = 2, k = 3)
  speed <- attr(r, "speed")
  expect_equal(unname(speed[c("R", "N")]),
               c(mean(web$rate_R[web$R > 0]), mean(web$rate_N[web$N > 0])))
  expect_true(is.na(speed[["P"]]))
  tip_edge <- match(1:3, tree$edge[, 2])
  expect_equal(r[tip_edge], unname(2 * (1 + 3 * speed[c("R", "N", "R")])))
  expect_equal(r[-tip_edge], 2)
  # unnamed levels follow the order of the tips; names are matched to the tips
  expect_equal(food_web_rates(web, tree, c("R", "N", "R")),
               food_web_rates(web, tree, c(c = "R", b = "N", a = "R")))
  # half of the branch of c
  r_half <- food_web_rates(web, tree, c("R", "N", "R"), k = 3, window = 1)
  expect_equal(r_half[tip_edge[3]], 1 + 3 * speed[["R"]] * 0.5)
  out <- trait_rates(tree, edge_rates = r)
  expect_equal(as.numeric(out$edge_rate), as.numeric(r))
})

test_that("food_web_rates() checks its input", {
  expect_error(food_web_rates(list(), tree, c("R", "N", "R")), "food_web_dynamics")
  expect_error(food_web_rates(web, list(), c("R", "N", "R")), "phylo")
  expect_error(food_web_rates(web, tree, c("R", "N")), "each of the 3 tips")
  expect_error(food_web_rates(web, tree, c("R", "N", "X")), "each of the 3 tips")
  expect_error(food_web_rates(web, tree, c(a = "R", b = "N", d = "R")), "differ")
  expect_error(food_web_rates(web, tree, c("R", "N", "P")), "level P is absent")
  expect_error(food_web_rates(web, tree, c("R", "N", "R"), k = -1), "`k`")
})
