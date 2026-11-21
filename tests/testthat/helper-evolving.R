# species j (index over all species ever) buds: new species g takes `amounts[[g]]` individuals
# (all, when NULL) in the sites `sites[[g]]`; a one-species wrapper of evolving_split_many()
evolving_split <- function(st, j, sites, names_new, amounts = NULL, jump_var = 0,
                           mutants = FALSE) {
  evolving_split_many(st, list(list(j = j, sites = sites, amounts = amounts, names = names_new,
                                    ids = st$k_all + seq_along(sites))), jump_var, mutants)
}
