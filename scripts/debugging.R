

df <- psych::bfi %>%
  dplyr::mutate(gender = factor(gender,
                         labels = c("Males", "Females")))

m <- "
A =~ A1 + A2 + A3 + A4 + A5
E =~ E1 + E2 + E3 + E4 + E5 + A5

"

fit1 <-lavaan::cfa(model = m,
            data = psych::bfi %>%
              dplyr::mutate(gender = factor(gender,
                                            labels = c("Males", "Females"))),
            #missing = "ml",
            meanstructure = T,
            ordered = F,
            parameterization = "theta")

fit2 <- lavaan::cfa(model = m,
            data = df,
            missing = "pairwise",
            # group = "gender",
            meanstructure = T,
            ordered = T)

fit3 <- lavaan::cfa(model = m,
            data = df,
            meanstructure = T,
            missing = "pairwise",
            group = "gender",
            ordered = T)

info <- model_info(fit3)
lv <- info$latent_variables
dat2 <- lavPredict_parallel(fit2)
dat3 <- lavPredict_parallel(fit3)
info2 <- model_info(fit2)
info3 <- model_info(fit3)

.augment_ordinal(fit2)
.augment_ordinal(fit3)

.augment_ordinal(dat = dat3,
                 fit = fit3,
                 info = info3) %>%
  view()

.augment_ordinal(fit2)

aug <- augment(fit1)

aug <- augment_ordinal(fit3)
aug %>%
  ggplot(aes(A, .yhat_A4)) +
  geom_point()

names(aug)

dat %>%
  select(all_of(lv)) %>%
  map(range, na.rm = TRUE)

augment_continuous(fit3, se = T) %>%
  View()
lavPredict(fit3, se = T)
lavPredict(fit3, se = "standard")

lavPredict_parallel_2(fit1,
                      se = TRUE,
                      R = 100)

lavPredict_parallel_2(fit2, se = T)

model_info(fit1)
lavaan::lavInspect(fit1, "estimator")

aug <- augment_ordinal_5(fit3, ci = TRUE,
                         ystar = TRUE, pr = TRUE,
                         resid = TRUE, yhat = TRUE)

names(df)
x <- df %>%
  select(.pr_X1__A3:.pr_X6__A3) %>% rowSums()

all(near(x, 1))
table(x)

model_info(fit3)


augment2(fit2) %>%
  ggplot(aes(A, .yhat_A1)) +
  geom_ribbon(aes(ymin = .yhat_lwr_A1, ymax = .yhat_upr_A1),
              fill = "grey")+
  geom_line() +
  facet_wrap(~group)


dat1 <- lavaan::lavInspect(fit1, "data")
dat2 <- lavaan::lavInspect(fit3, "data")

group_var1    <- tryCatch(lavaan::lavInspect(fit1, "group"),       error = function(e) NULL)
group_labels1 <- tryCatch(lavaan::lavInspect(fit1, "group.label"), error = function(e) NULL)

group_var2    <- tryCatch(lavaan::lavInspect(fit2, "group"),       error = function(e) NULL)
group_labels2 <- tryCatch(lavaan::lavInspect(fit2, "group.label"), error = function(e) NULL)


dat2 <- dat2 %>%
  map(as_tibble) %>%
  bind_rows(.id = group_var2)

# Comments are in English.
library(future)
library(future.apply)
library(furrr)

# 1) Choose a plan (multisession works cross-platform)
workers <- parallel::detectCores() - 1
plan(multisession, workers = workers)

unique_valid <- function(x){
  out <- unique(x)
  out[!is.na(out)]
}

create_dummy <- function(data) {
  values <- data %>%
    map(unique_valid)
  longest <- values %>%
    map_int(length) %>%
    max()
  values %>%
    map(rep_len,
        length.out = longest) %>%
    bind_cols()
}

# 2) Split data into chunks (by rows) — or by groups for multi-group models
dat_original <- lavInspect(fit1, "data") %>%
  as_tibble()

dat_unique <- dat_original %>%
  distinct()

ov_ord <- lavNames(fit1, "ov.ord")
n_rows <- nrow(dat)
n <- n_rows / workers


idx_list <- split(seq_len(n_rows),
                  ceiling(seq_along(seq_len(n_rows))/n))

chunked <- vector("list", length = length(idx_list))

for (i in seq_along(idx_list)) {
  chunked[[i]] <- dat_unique[idx_list[[i]], , drop = FALSE]
}


dummy <- dat_unique %>%
  create_dummy()

chunked <- chunked %>%
  future_map(~bind_rows(dummy, .x))

out <- chunked %>%
  future_map(~lavPredict(fit1, newdata = .x,
                         append.data = TRUE,
                         assemble = TRUE))

dummy_rows <- nrow(dummy)

for (i in seq_along(out)) {
  out[[i]] <- out[[i]][-c(seq(dummy_rows)) , ]
}

out <- do.call(rbind, out) %>%
  as_tibble()

dat_original %>%
  left_join(out)

list(out)

# 2) Split data into chunks (by rows) — or by groups for multi-group models
group_var  <- tryCatch(lavaan::lavInspect(fit3, "group"), error = function(e) NULL)
ov_ord <- lavNames(fit3, "ov.ord")

dat_original <- lavInspect(fit3, "data") %>%
  map(as_tibble)

group_labels  <- tryCatch(lavaan::lavInspect(fit3, "group.label"), error = function(e) NULL)

dat_original <- dat_original %>%
  bind_rows(.id = group_var)

dat_unique <- dat_original %>%
  distinct()

dummy <- dat %>%
  create_dummy()

n_rows <- nrow(dat_unique)
n <- n_rows / workers


idx_list <- split(seq_len(n_rows),
                  ceiling(seq_along(seq_len(n_rows))/n))

chunked <- vector("list", length = length(idx_list))

for (i in seq_along(idx_list)) {
  chunked[[i]] <- dat_unique[idx_list[[i]], , drop = FALSE]
}

chunked

chunked <- chunked %>%
  future_map(~bind_rows(dummy, .x))

out <- chunked %>%
  future_map(~lavPredict(fit3,
                         newdata = .x,
                         append.data = TRUE,
                         assemble = TRUE,
                         drop.list.single.group = FALSE))

dummy_rows <- nrow(dummy)

for (i in seq_along(out)) {
  out[[i]] <- out[[i]][-c(seq(dummy_rows)) , ]
}

out <- do.call(rbind, out) %>%
  as_tibble()

out <- dat_original %>%
  left_join(out)

group_col <- out[, group_var]

out[, group_var] <- NULL

out <- split(out, group_col)




out[group_labels]





out <- out %>%
  future_map(as_tibble) %>%
  future_map(~slice(.x, -(1:dummy_rows))) %>%
  bind_rows()

group_col <- out[[group_var]]

out[[group_var]] <- NULL

split(out, group_col)


out <- out %>%
  mutate(
    across(everything(), as.double)
  )

lavaan::lavPredict(
  fit, transform = FALSE, append.data = TRUE,
  assemble = FALSE, drop.list.single.group = FALSE
)

dummy

n_chunks <- 10

ov_ord <- lavNames(fit1, "ov.ord")

n_rows <- nrow(dat)

n <- n_rows / n_chunks


split_data <- function(data, n_chunks = 10) {
  n_rows <- nrow(data)
  n <- n_rows / n_chunks


  idx_list <- split(seq_len(n_rows),
                    ceiling(seq_along(seq_len(n_rows))/n))

  chunked <- vector("list", length = length(idx_list))

  for (i in seq_along(idx_list)) {
    chunked[[i]] <- data[idx_list[[i]], , drop = FALSE]
  }

  return(chunked)

}


dat <- dat %>%
  map(split_data)



chunked <- vector("list", length = length(group_labels)) %>%
  set_names(group_labels)

for (group in group_labels) {
  chunks <- dat[[group]]
  padding <- dummy[[group]]

  for (i in seq_along(chunks)) {
    chunked[[group]][[i]] <- bind_rows(padding, chunks[[i]])
    chunked[[group]][[i]][[group_var]] <- group
  }

}

map(chunked,
    ~map(
      .x,
      ~lavPredict(fit3,
                  newdata = .x,
                  append.data = TRUE)
    ))

chunked %>%
  map(
    ~map(~lavPredict(fit3, newdata = .x,
                    append.data = TRUE))
  )

out <- chunked %>%
  map(~lavPredict(fit3, newdata = .x,
                         append.data = TRUE))
out

d
dat[[2]]

dat %>%
  map(length)



dat <- dat %>%
  bind_rows(.id = group_var ) %>%
  mutate(chunk = rep_len(1:10, length.out = nrow(.))) %>%
  group_by(across(all_of(c(group_var, "chunk")))) %>%
  summarise(
    # keep all columns of the current group (incl. grouping vars) in a list-col
    data = list(cur_data_all()),
    .groups = "drop"
  ) %>%
  left_join(dummy, by = group_var)



dat <- dat %>%
  mutate(
    data = map2(dummy, data, bind_rows),
    n_drop = map_int(dummy, nrow)
  ) %>%
  select(-dummy)

dat <- dat %>%
  mutate(
    data = map(data,
               ~select(.x, -all_of(c("chunk", group_var))))
  )

dat <- dat %>%
  unnest(data) %>%
  nest(data = all_of(c(ov_ord, group_var)))

dat <- dat %>%
  mutate(
    fs = future_map(
      data,
      ~lavPredict(fit3,
                  newdata = .x,
                  append.data = TRUE)
    )
  )

convert <-

dat %>%
  mutate(

  )


dat %>%
  mutate(
    fs = map(
      fs,
      ~.x %>%
        as_tibble() %>%
        bind_rows(.id = group_var))
  )

dat <- dat %>%
  mutate(chunk = rep_len(1:10, length.out = nrow(.)))


dummy <- dat %>%
  map(
    ~(.x %>%
        select(all_of(ov_ord)) %>%
        future_map(unique_valid) %>%
        bind_cols() %>%
        tidyr::fill(everything(), .direction = "downup"))

  )

group_var  <- tryCatch(lavaan::lavInspect(fit3, "group"),       error = function(e) NULL)

make_chunks <- function(data) {
  idx_list <- split(seq_len(nrow(data)), ceiling(seq_along(seq_len(nrow(data)))/100))

  print(idx_list)

  chunked <- vector("list", length = length(idx_list))

  for (i in seq_along(idx_list)) {
    chunked[[i]] <- data[idx_list[[i]], , drop = FALSE]
  }

  return(chunked)
}

chunked <- dat %>%
  map(make_chunks)


for (group in names(chunked)) {
  map(
    chunked[[]]
  )

}

nms <- names(chunked)

seq_along(chunked)

names(chunked)

out <- vector("list", length = length(chunked))

for (group in names(chunked)) {
  for (chunk in chunked[[group]]) {
    print(str_c(group, chunk))
  }
}

chunk
seq_along(chunked)

str()

for (i in seq_along(chunked)) {
  for (j in vector) {

  }
}

chunked %>%
  map(
    ~map2(

    )
  )

map2(
  dummy,
  chunked,
  ~map(
    ~bind_rows(.x, .y)
  )
)



str(dat)

dat <- dat %>%
  bind_rows(.id = group_var)

idx_list <- split(seq_len(nrow(dat)), ceiling(seq_along(seq_len(nrow(dat)))/100))

chunked <- vector("list", length = length(idx_list))

for (i in seq_along(idx_list)) {
  chunked[[i]] <- dat[idx_list[[i]], , drop = FALSE]
}

str(chunked)

ov_ord <- lavNames(fit3, "ov.ord")

unique_valid <- function(x){
  out <- unique(x)
  out[!is.na(out)]
}

dummy <- dat %>%
  select(all_of(ov_ord)) %>%
  future_map(unique_valid) %>%
  bind_cols() %>%
  tidyr::fill(everything(), .direction = "downup")


chunked <- chunked %>%
  future_map(~bind_rows(dummy, .x))


lavPredict(fit1, chunk_)

parallel_lavPredict(fit1, chunk_size = 200)

# 1) Choose a plan (multisession works cross-platform)
plan(multisession, workers = parallel::detectCores() - 1)

# 2) Split data into chunks (by rows) — or by groups for multi-group models
dat <- lavInspect(fit1, "data") %>%
  as_tibble()

ov_ord <- lavNames(fit1, "ov.ord")

dat <- dat %>%
  mutate(
    across(all_of(ov_ord), as.ordered)
  )

idx_list <- split(seq_len(nrow(dat)), ceiling(seq_along(seq_len(nrow(dat)))/100))

chunked <- vector("list", length = length(idx_list))

for (i in seq_along(idx_list)) {
  chunked[[i]] <- dat[idx_list[[i]], , drop = FALSE]
}

str(chunked)

plan("multisession", workers = 10)

out <- chunked %>%
  future_map(~lavPredict(fit1, newdata = .x,
                         append.data = TRUE))

n_rows <- nrow(dummy)

out <- out %>%
  future_map(as_tibble) %>%
  future_map(~slice(.x, -(1:n_rows))) %>%
  bind_rows()


parallel_lavPredict(fit1)

plot_cfa(fit1)
