data("variants")

context("VERSO")

data("inference")
test_that("VERSO produces correct output", {
    expect_equal(names(inference),c("B","C","phylogenetic_tree","corrected_genotypes","genotypes_prevalence","genotypes_summary","log_likelihood","error_rates"))
})
