




### Criteria #3 Monotony just with the log2(OR)

or1 <- OR_df[,2]
or2 <- OR_df[,3]
or3 <- OR_df[,4]

### Verify if terms are decreasing monotone 
decreasing_monotony <- OR_df[(or1 > or2 & or2 > or3),]

### Verify if terms have not trend 
equal_trend <- OR_df[(or1 == or2 & or2 == or3),]

### Verify if terms are increasing monotone 
increasing_monotony <- OR_df[(or1 < or2 & or2 < or3),]

dir.create(path = file.path(output_dir,"9_Trend_criteria"))
write.table(decreasing_monotony,file.path(output_dir,"9_Trend_criteria","decreasing_monotony_3rd_criteria.txt"))
write.table(increasing_monotony,file.path(output_dir,"9_Trend_criteria","increasing_monotony_3rd_criteria.txt"))



# 1 load dfs and merge them



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
df$Pattern_C2[which(df$OR_GD>df$upper_CI_G & (((df$lower_CI_G > df$lower_CI_GL) & (df$upper_CI_G < df$upper_CI_GL))|
                                                ((df$lower_CI_G < df$lower_CI_GL) & (df$upper_CI_G > df$upper_CI_GL))))]=2
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
summary(factor(df$Pattern_C2))









