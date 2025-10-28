

df <- psych::bfi %>%
  dplyr::mutate(gender = factor(gender,
                         labels = c("Males", "Females")))

m <- "
A =~ A1 + A2 + A3 + A4 + A5
E =~ E1 + E2 + E3 + E4 + E5
"

ord <- c("E1", "E2", "E3", "E4", "E5")

fit1 <-lavaan::cfa(model = m,
            data = psych::bfi %>%
              dplyr::mutate(gender = factor(gender,
                                            labels = c("Males", "Females"))),
            #missing = "ml",
            meanstructure = T,
            ordered = T,
            parameterization = "theta")

fit1 <- lavaan::cfa(model = m,
            data = df,
            missing = "pairwise",
            # group = "gender",
            meanstructure = T,
            ordered = F)

fit2 <- lavaan::cfa(model = m,
                    data = df,
                    meanstructure = T,
                    missing = "pairwise",
                    group = "gender",
                    ordered = ord)
resid_cor(fit1)
resid_cor(fit2)

lavaan::lavResiduals(fit1)

# single group
plot_cfa(fit1)
resid_corrplot(fit1)

resid_qq(fit1)

lavPredict_parallel(fit1)
lavPredict_parallel(fit3)

model_info(fit3)


lavaan::lavInspect(fit3, "case.idx")
.case_ids <- function(fit) {
  lavaan::lavInspect(fit, "case.idx")
}

lavaan::lavPredict(fit2, append.data = TRUE, drop.list.single.group = FALSE)
lavaan::lavPredict(fit3, append.data = TRUE, drop.list.single.group = FALSE)

hopper_plot(fit2)
hopper_plot(fit3)
.prepare_continuous(fit2)
.prepare_ordinal(fit2)

.prepare_continuous(fit3)
.prepare_ordinal(fit3)


lavPredict_parallel(fit2)
lavPredict_parallel(fit3)

sg1 <- item_data(fit1) # single group model
sg2 <- item_data(fit2) # single group model


item_plot(sg1, latent = "E",
          sort = "r2",
          jitter_seed = 333)
item_plot(sg1, latent = "A",
          sort = "r2",
          jitter_seed = 333)
item_plot(sg2, latent = "E",
          sort = "r2",
          jitter_seed = 333)
item_plot(sg2, latent = "A",
          sort = "r2",
          jitter_seed = 333)

p2 <- item_plot(sg2, latent = "E",
                sort = "r2",
                jitter_seed = 333)
print(p2)

p2 <- item_plot(mg,
                latent = "A",
                facet = "grid",
                metrics_pad = .1,
                alpha_points = .1,
                sort = "r2",
                jitter_sd = 0.1)
print(p2)

non_pen <- c("r2", "rmse", "mae")
pen <- c("r2_pen", "rmse_pen", "mae_pen")

facet_names <- c(
  r2   = "1 - italic(R)^2",  # R kurzívou a 2 jako horní index
  rmse = "\"RMSE\"",
  mae  = "\"MAE\""
) %>%
  ggplot2::as_labeller(ggplot2::label_parsed)

sg1$metrics %>%
  mutate(
    across(
      c(r2, r2_pen),
      ~1 - .x
    )) %>%
  pivot_longer(cols = all_of(non_pen),
               names_to = "metric") %>%
  mutate(
    metric = factor(metric, levels = non_pen),
    item = reorder_within(item, -value, within = metric)) %>%
  ggplot(aes(value, item)) +
  geom_point(size = 2) +
  geom_line(mapping = aes(group = 1)) +
  facet_wrap(~metric, scales = "free",
             labeller = facet_names) +
  ggplot2::expand_limits(x = 0) +
  ggplot2::scale_y_discrete(labels = function(x) sub(paste0("__", ".*$"), "", x)) +
  labs(x = "Value", y = "Item")

p <- position_dodge(width = .2, orientation = "y")

sg2$metrics %>%
  mutate(
    across(
      c(r2, r2_pen),
      ~1 - .x
    )) %>%
  pivot_longer(cols = all_of(non_pen),
               names_to = "metric") %>%
  mutate(
    metric = factor(metric, levels = non_pen),
    item = reorder_within(item, -value, within = metric)) %>%
  ggplot(aes(value, item, color = .group)) +
  geom_point(size = 2, position = p) +
  geom_line(mapping = aes(group = .group),
            orientation = "y",
            position = p) +
  facet_wrap(~metric, scales = "free",
             labeller = facet_names) +
  ggplot2::expand_limits(x = 0) +
  ggplot2::scale_y_discrete(labels = function(x) sub(paste0("__", ".*$"), "", x)) +
  labs(x = "Value", y = "Item", color = "Group") +
  theme(legend.position = "top")


# Single group model
orig <- sg$original_data
metrics <- sg$metrics
new <- sg$new_data

metrics <- metrics %>%
  mutate(
    across(
      c(r2, rmse, mae, r2_pen, rmse_pen, mae_pen),
      ~.x %>%
        round(digits = 3) %>%
        format(digits = 3, nsmall = 3)
    )
  )

item <- metrics %>%
  filter(item == "A2")


r2   <- item$r2;   r2_pen   <- item$r2_pen
rmse <- item$rmse; rmse_pen <- item$rmse_pen
mae  <- item$mae;  mae_pen  <- item$mae_pen

cap <- bquote(
  italic(R)^2 == .(r2) * "," ~
    plain("RMSE") == .(rmse) * "," ~
    plain("MAE") == .(mae) * ","
)

set.seed(123)

n <- 100              # počet respondentů
k <- 3                # počet položek
lambda <- c(.4, .5, .6)

F1 <- rnorm(n)        # latentní faktor
E  <- matrix(rnorm(n * k), n, k)  # náhodné chyby
Y  <- sweep(E, 2, lambda * F1, "+")  # Y_ij = λ_j * F_i + ε_ij

colnames(Y) <- paste0("y", 1:k)
Y <- as.data.frame(Y)
head(Y)

congen <- psych::sim.congeneric(
  loads = seq(0.4, 0.6, by = .1),
  N = 1000,
  short = FALSE,
  categorical = TRUE,
  cuts = 0
) %>%
  .$observed %>%
  dplyr::as_tibble() %>%
  rlang::set_names("a1", "a2", "a3")




mod <- "
F =~ V1 + V2 + V3 + V4 + V5 + V6
"

fit5 <-lavaan::cfa(model = mod,
                   data = congen,
                   ordered = T)



idata <-  item_data(fit5)


item_plot(idata, latent = "F",
          sort = "r2",
          jitter_seed = 333)

.prepare_continuous(fit5, other_latents = "mean")

.prepare_ordinal(fit5)

psych::con


load("data/ipip.RData")
df

table
names(table)
