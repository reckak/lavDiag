

df <- psych::bfi %>%
  dplyr::mutate(gender = factor(gender,
                         labels = c("Males", "Females")))

m <- "
A =~ A1 + A2 + A3 + A4 + A5
E =~ E1 + E2 + E3 + E4 + E5 + A5

"

fit1 <-lavaan::cfa(model = m,
            data = df,
            missing = "ml",
            meanstructure = F)

fit2 <- lavaan::cfa(model = m,
            data = df,
            missing = "ml",
            group = "gender",
            ordered = F)

fit3 <- lavaan::cfa(model = m,
            data = df,
            missing = "pairwise",
            group = "gender",
            ordered = TRUE)

resid_cor(fit1)

is_not_lavaan_fit(1L)

hopper_plot(fit3)
