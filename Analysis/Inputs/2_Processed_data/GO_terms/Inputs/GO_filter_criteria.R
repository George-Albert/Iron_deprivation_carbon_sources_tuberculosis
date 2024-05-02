
df_up=read.table("Inputs/df_up_GO_enrichment.txt")

df_down=read.table("Inputs/df_down_GO_enrichment.txt")

list_up <- readRDS("Inputs,lista_up_GO_enrichment.rds")
list_down <- readRDS("Inputs,lista_down_GO_enrichment.rds")

### Save the max and min values distinct from Inf
max_values <- apply(df[,2:ncol(df)], 2, function(x) max(x[is.finite(x)], na.rm = TRUE))
min_values <- apply(df[,2:ncol(df)], 2, function(x) min(x[is.finite(x)], na.rm = TRUE))

for (col in colnames(df[,2:ncol(df)])) {
  df[is.infinite(df[, col]), col] <- max_values[col] * 10
}


## 2. add flag columns

df$Pattern_C1=0

df$Pattern_C2=0

df$Pattern_C3=0

### Pattern #1
### Set 1: bajar, bajar:
df$Pattern_C1[which(df$lower_CI_GD>df$upper_CI_G & df$lower_CI_G > df$upper_CI_GL)]=1
## 1
## 2
## 3
## set 2: bajar, stat
df$Pattern_C1[which(df$lower_CI_GD>df$upper_CI_G & (((df$lower_CI_G > df$lower_CI_GL) & (df$lower_CI_G < df$upper_CI_GL))|
                                                      ((df$upper_CI_G > df$lower_CI_GL) & (df$upper_CI_G < df$upper_CI_GL))))]=2
## 1
## 2
## 3
### Set 1: subir, subir:
df$Pattern_C1[which(df$upper_CI_GD < df$lower_CI_G & df$upper_CI_GD < df$lower_CI_GL)]=3
## 1
## 2
## 3

### Pattern #2
### Set 1: bajar, bajar
df$Pattern_C2[which(df$OR_GD > df$upper_CI_G & df$OR_G > df$upper_CI_GL)]=1
### Set 2: bajar, stat
df$Pattern_C2[which(df$OR_GD>df$upper_CI_G & (((df$lower_CI_G > df$lower_CI_GL) & (df$lower_CI_G < df$upper_CI_GL))|
                                                ((df$upper_CI_G > df$lower_CI_GL) & (df$upper_CI_G < df$upper_CI_GL))))]=2
### Set 3: bajar, subir
df$Pattern_C2[which(df$OR_GD > df$upper_CI_G & df$OR_G < df$upper_CI_GL)]=3

### Set 4: stat, bajar
### Set 5: stat, stat
### Set 6: stat, subir

### Pattern #3

### Set 1: bajar, bajar
df$Pattern_C3[which(df$OR_GD > df$OR_G & df$OR_G > df$OR_GL)]=1
### Set 2: bajar, stat
df$Pattern_C3[which(df$OR_GD > df$OR_G & df$OR_G == df$OR_GL)]=2
### Set 3: bajar, subir
df$Pattern_C3[which(df$OR_GD > df$OR_G & df$OR_G < df$OR_GL)]=3


### Set 4: stat,bajar
df$Pattern_C3[which(df$OR_GD == df$OR_G & df$OR_G > df$OR_GL)]=4
### Set 5: stat,stat
df$Pattern_C3[which(df$OR_GD == df$OR_G & df$OR_G == df$OR_GL)]=5
### Set 6: stat,subir
df$Pattern_C3[which(df$OR_GD == df$OR_G & df$OR_G < df$OR_GL)]=6

### Set 7: subir, bajar
df$Pattern_C3[which(df$OR_GD < df$OR_G & df$OR_G > df$OR_GL)]=7
### Set 8: subir,stat
df$Pattern_C3[which(df$OR_GD < df$OR_G & df$OR_G == df$OR_GL)]=8
### Set 9: subir, subir
df$Pattern_C3[which(df$OR_GD < df$OR_G & df$OR_G < df$OR_GL)]=9


# Count Flags
summary(factor(df$Pattern_C3))

write.table(df,file=file.path(input_folder,"df_up_to_now.txt"))
