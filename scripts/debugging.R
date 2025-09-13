

df <- psych::bfi %>%
  dplyr::mutate(gender = factor(gender,
                         labels = c("Males", "Females")))

m <- "
A =~ A1 + A2 + A3 + A4 + A5
E =~ E1 + E2 + E3 + E4 + E5 + A5

"

fit1 <-lavaan::cfa(model = m,
            data = df,
            # missing = "ml",
            meanstructure = T)

fit2 <- lavaan::cfa(model = m,
            data = df,
            missing = "ml",
            group = "gender",
            ordered = F)

fit3 <- lavaan::cfa(model = m,
            data = df,
            meanstructure = T,
            missing = "pairwise",
            group = "gender",
            ordered = F)

resid_cor(fit1)

lavDiag::.is_not_lavaan_fit(1L)

hopper_plot(fit1)

resid_qq(fit2, n = 5)

resid_corrplot(fit2)

model_info(fit2)

parameter_estimates(fit2)

plot_cfa(fit1)

info <- model_info(fit1)

fs_and_ov <- lavaan::lavPredict(fit1,
                                transform = F,
                                append.data = T,
                                assemble = T,
                                drop.list.single.group = TRUE)

eta_hat <- fs_and_ov[, info$latent_variables]

params <- lavInspect(fit1, "est")
lambda <- params$lambda
nu <- params$nu

y_hat <- sweep(eta_hat %*% t(lambda), 2, nu, FUN = "+")
nms <- colnames(y_hat)
nms_yhat <- str_c(".yhat_", nms)

colnames(y_hat) <- nms_yhat
y_hat

out <- fs_and_ov %>%
  dplyr::bind_cols(y_hat)

ov <- info$observed_variables

for (nm in ov) {
  observed <- nm
  yhat <- str_c(".yhat_", nm)
  resid <- str_c(".resid_", nm)
  out[[resid]] <- out[[observed]] - out[[yhat]]
}

fs_and_ov <- lavaan::lavPredict(fit3,
                                transform = F,
                                append.data = T,
                                assemble = F,
                                drop.list.single.group = FALSE)
fs_and_ov

params <- lavInspect(fit3, "est")

info <- model_info(fit3)

eta_hat <- fs_and_ov %>%
  map(~.x[, info$latent_variables])

lambda <- params %>%
  map("lambda")

nu <- params %>%
  map("nu")

y_hat <- vector("list", info$n_groups) %>%
  setNames(info$group_labels)

for (i in seq_len(info$n_groups)) {
  y_hat[[i]] <- sweep(eta_hat[[i]] %*% t(lambda[[i]]), 2, nu[[i]], FUN = "+")
}


for (i in seq_len(info$n_groups)) {
  nms <- colnames(y_hat[[i]])
  new_names <- str_c(".yhat_", nms)
  colnames(y_hat[[i]]) <- new_names
}


y_hat <- y_hat %>%
  map(as_tibble) %>%
  bind_rows()

fs_and_ov <- fs_and_ov %>%
  map(as_tibble) %>%
  bind_rows(.id = "group")

fs_and_ov <- fs_and_ov %>%
  mutate(
    across(-group, as.double)
  )

out <- fs_and_ov %>%
  dplyr::bind_cols(y_hat)

ov <- info$observed_variables

for (nm in ov) {
  observed <- nm
  yhat <- str_c(".yhat_", nm)
  resid <- str_c(".resid_", nm)
  out[[resid]] <- out[[observed]] - out[[yhat]]
}

augment(fit3)

lavaan::inspect(fit3, "lambda")
