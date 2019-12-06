## Análise em Nível de Transcrito Kallisto
## Utilizando Octuplicatas
## Organização do data frame

# Nomes de populações
ZIKA <- 'ZIKA'
CHIKV <- 'CHIKV'
CHIKV_REC <- 'CHIKV_REC'
GBS <- 'GBS'
GBS_REC <- 'GBS_REC'
CONTROL <- 'CONTROL'


# Vetor com nomes populações para fazer coluna pop:
# Obs.: Ao fazer o vetor, melhor deixar organizado combinando amostra + nomearquivo
pop <- c(rep(GBS, 8), rep(CONTROL, 16), rep(ZIKA, 8), rep(GBS, 16), 
         rep(GBS_REC, 8), rep(CHIKV_REC, 8), rep(CONTROL, 8), rep(GBS, 8), 
         rep(CONTROL, 8), rep(CHIKV, 8), rep(ZIKA, 8), rep(GBS, 8), 
         rep(GBS_REC, 24), rep(CHIKV, 16), rep(CONTROL, 16), rep(CHIKV_REC, 8), 
         rep(ZIKA, 8), rep(GBS_REC, 16), rep(GBS, 8), rep(CHIKV, 8), 
         rep(CHIKV_REC, 32), rep(ZIKA, 8), rep(GBS_REC, 24))


pop          # ZIKA, CHIKV, CHIKV_REC, GBS, GBS_REC, CONTROL
length(pop)  # 280


# Nome de centro de pesquisa para fazer coluna center
center <- rep('IMT-UFRN', 280)
center
length(center) # 280 == pop


# Nomes de amostras analisadas para fazer coluna run
run <- list.files(dir)
run
mode(run)
length(run)  # 280 == pop


condition <- c(rep('gbs', 8), rep('control', 16), rep('zika', 8), rep('gbs', 16), 
               rep('gbs_rec', 8), rep('chikv_rec', 8), rep('control', 8), rep('gbs', 8), 
               rep('control', 8),rep('chikv', 8), rep('zika', 8), rep('gbs', 8), 
               rep('gbs_rec', 24), rep(CHIKV, 16), rep('control', 16), rep('chikv_rec', 8), 
               rep('zika', 8), rep('gbs_rec', 16), rep('gbs', 8), rep(CHIKV, 8), 
               rep('chikv_rec', 32), rep('zika', 8), rep('gbs_rec', 24))

length(condition)

# Replicatas
replicates <- c('rep01', 'rep02', 'rep03 ', 'rep04',
                'rep05', 'rep06', 'rep07 ', 'rep08')


# Unir cada vetor formando colunas de um data frame:
samples_info <- data.frame(pop = pop, 
                      center = center, 
                      run = run,
                      condition = condition,
                      replicate = rep(replicates, 35))   


head(samples_info, 10)
str(samples_info)
names(samples_info)
length(samples_info$replicate)  # 280


