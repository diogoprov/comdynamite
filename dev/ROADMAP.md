# comdynamite: estado e próximos passos

Atualizado em 2026-11-06. Versão atual: 0.6.0 (DOI de conceito `10.5281/zenodo.23993533`).
Tags: `v0.1.0`, `v0.2.0`, `v0.5.0`, `v0.5.1`, com o que faltava da v0.5.0, `v0.5.2` e `v0.6.0`. `v0.3.0` e `v0.4.0`
não têm tag porque seguem em aberto: a ordem mudou em 2026-11-06, quando os métodos da
v0.5.0 passaram à frente, e a próxima versão depois deles será a do marco que fechar.
As versões abaixo são uma sugestão: cada marco tem um critério de conclusão, e a
ordem pode mudar. Tag no formato `vX.Y.Z`, com release no GitHub.

## Visão

O comdynamite deixa de ser só um simulador e passa a ser um hub de métodos que ligam
traits, filogenias, biogeografia e metacomunidades. A simulação continua no centro:
cada método do pacote pode ser testado com dados em que a verdade é conhecida.
Métodos publicados que nunca chegaram ao R, como os de Leibold et al. (2010) e
Peres-Neto et al. (2012), são pouco usados por isso e têm prioridade.

Para isso funcionar, o pacote precisa de uma estrutura de dados comum, dependências
leves (pacotes grandes em `Suggests`) e convenção de nomes (`simulate_*` para gerar,
e um prefixo a definir para analisar). O título e o `Description` mudaram na 0.6.0.

## Feito

| Item | O que entrou |
|---|---|
| Modelo | Thompson et al. (2020) com traits e competição filogenéticos; paisagens implícitas e em grade; permutações do toro; ambiente variável no tempo |
| Montagem por regras | `assemble_community()`: similaridade limitante, filtragem e neutra |
| Traits | `traits_ou()` (um ou vários ótimos) e `paint_regimes()`, no lugar de `traits_bm()` |
| Identidade | hex sticker (opção F1) e `data-raw/hex_sticker.R` |
| Site | pkgdown em https://provetelab.org/comdynamite/, com build automático no `main` |
| Simulador estático | `simulate_static()`, na versão de desenvolvimento |
| Vignettes | `comdynamite`, `phylogenetic-structure` e `evaluating-methods` |
| Documentação | `\value` e exemplos em todas as funções exportadas; `Description` atualizada |

## v0.1.0: primeira versão citável (concluída em 2026-11-05)

CI em macOS, Windows e Linux; 100% de cobertura de linhas; `inst/CITATION` e
`CITATION.cff`; release no GitHub com DOI do Zenodo.

- [x] Submissão ao r-universe: `comdynamite` incluído em `packages.json`, no repositório
      `diogoprov/diogoprov.r-universe.dev`.
- [ ] Confirmar o build em https://diogoprov.r-universe.dev/comdynamite. Depois, pôr o
      badge e a instalação pelo r-universe no README, como no treeshapesim.

## v0.2.0: validação e geradores estáticos

Objetivo: mostrar que a implementação reproduz o modelo publicado e ter uma linha de
base sem dinâmica.

- [x] `simulate_static()`: resposta gaussiana de Dray & Legendre (2008), usada também
      por Peres-Neto et al. (2012). Nessas fontes, `h` vem de U(0,5; 1), a abundância é
      o valor esperado e a incidência é binomial com probabilidade `y`, opções de
      `sampling`. As contagens de Poisson e `h` em U(0, 30), padrão no pcps-sim, vêm de
      Peres-Neto & Kembel (2015, Appendix A) e Duarte et al. (2016, Appendix S2). Duarte
      et al. atribuem o Poisson a 2012, mas ele não está no texto principal desse artigo.
- [x] Dois gradientes em `simulate_static()`: produto das respostas ao ambiente e ao
      espaço, como em `ScriptDemo.m`.
- [x] Gerador logístico do Appendix C (Ives & Helmus 2011), com os sete cenários de
      [E] e [S]: `simulate_logistic()` e `logistic_gradients()`.
- [x] Comparar com o `mcomsimr` nos mesmos parâmetros, levando em conta as diferenças
      de convenção (`sigma`, `alpha[j, k]`). Em `data-raw/validation_mcomsimr.R`: 16
      cenários (dispersão, competição, paisagem), 25 sítios, 12 espécies, 300 passos,
      ambiente constante e 200 réplicas por modelo. Nenhuma das 80 comparações de
      distribuição foi significativa (KS, FDR 5%), e as diferenças nas médias de alfa,
      gama e abundância total ficaram abaixo de 1%. Com extirpação de 0,02, nenhuma das 80
      foi significativa, mas alfa ficou 8% abaixo na grade com dispersão 0,001: nossos
      imigrantes são Poisson, como no código Julia dos autores, e o mcomsimr conserva os
      dispersores (alfa 1,75, contra 1,85 com multinomial e 1,85 no mcomsimr). Falta um
      ambiente que mude no tempo.
- [x] Reproduzir os padrões de Thompson et al. (2020), em
      `data-raw/validation_thompson.R` e `validation_thompson_cc.R`: 100 manchas, 50
      espécies, 2200 passos, 30 paisagens, 15 taxas de dispersão, 6 das 13 larguras de
      nicho, 4 tipos de competição. Os scripts seguem o código dos autores
      (Meta_com_framework v1.0, Zenodo 10.5281/zenodo.3708350), que difere do artigo em
      pontos que mudam o resultado: beta é gama/alfa, e não o gama − alfa da Tabela 1;
      manchas e tempos vazios não entram em alfa; o ambiente é convertido em postos,
      portanto uniforme em [0, 1]; dispersão começa em 0,00001, não 0,0001; a
      competição-colonização reduz os imigrantes à metade. Com isso, em dispersão 1e-5,
      stabilising dá alfa 4,6, gama 42 e beta 9,1 (σ = 0,5) e 6,6, 46 e 7,0 (σ = 10),
      próximos das Figuras 3c e 3d lidas a olho. H1, H2, H3, H5 e H6 aparecem.
- [x] H4 resolvido. A diferença anterior vinha do ambiente, não da definição de beta,
      que já era gama/alfa. A parte da competição aparece (beta espacial em 1e-5, σ = 0,5:
      equal 38,8, mixed 21,5, stabilising 9,1). A parte do nicho depende do cenário, como
      na Figura 4 deles: em stabilising, beta tem pico em σ = 0,215 (14,4 contra 7,0 em
      σ = 10); em mixed o pico é fraco, em σ = 0,5; em equal e dispersão baixa, beta cresce
      com o nicho (16,5, 38,8 e 43,1 em σ = 0,215, 0,5 e 10), por deriva. Nichos estreitos
      só elevam beta em equal a partir de dispersão perto de 0,02. O texto de H4 generaliza
      mais do que a figura mostra. A comparação com as figuras é visual, porque os valores
      publicados não estão em tabela.
- [x] `extirpation` em `simulate_metacommunity()` e `mc_dynamics()`. O cenário
      competição-colonização lembra as Figuras 3g e 3h: com σ = 0,5, gama sobe de 3,6 para
      cerca de 6 perto de dispersão 5e-4; com σ = 10, colapsa em uma espécie perto de 0,02.
- [x] Testes de propriedades com resultado analítico, em `test-analytic.R`: equilíbrios
      com uma e várias espécies, exclusão, fonte e sumidouro, fluxo pelo kernel.
- [x] Artigo "Validating the model of Thompson et al. (2020)" no site, com as divergências
      entre artigo e código, também resumidas no help de `simulate_metacommunity()`.

**Concluído quando:** houver um script em `data-raw/` que reproduza a comparação.

## Métodos de análise

- [x] `pcc_test()` e `phylo_eigenvectors()`: PCC de Peres-Neto & Kembel (2015), vindo do
      pcps-sim, com permutações restritas opcionais. Em `data-raw/validation_pcc.R`, as
      taxas de neither e environment batem com a Tabela 1 (erro tipo I do teste de
      espécies 0,049 e 0,040 contra 0,054 e 0,052). O poder (0,87 e 0,85 contra 0,945 e
      0,947) e o teste de espécies em phylogeny (0,20 e 0,19 contra 0,27 e 0,28) ficam
      abaixo, como esperado sem a árvore deles.
- [x] `phylo_beta_test()` e `phylo_beta()`: framework de beta diversidade filogenética do
      mesmo artigo (db-RDA sobre MPD entre comunidades), que não existia no pcps-sim. F e
      R2 iguais aos do `vegan::adonis2()`, e métrica igual a `picante::comdist()`. Em
      `data-raw/validation_phylo_beta.R`, com 1000 réplicas e o PCC nos mesmos dados, as
      taxas batem com a Tabela 1: neither .005 e .010 (.007, .008), environment .050 e
      .062 (.055, .050), phylogeny .033 e .056 (.035, .039), teste de espécies em
      phylogeny .16 e .19 (.167, .177). Poder .893 e .908 (.937, .939), R2 ajustado .083
      e .060 (.091, .064). O PCC teve poder de .888 e .856 nos mesmos dados, ao contrário
      do artigo, onde fica à frente (.945, .947), e um erro de .058 onde a outra rodada
      deu .068, o que mostra o ruído desses valores.
- [ ] Heatmaps do artigo para visualizar os eixos de `scores`.
- [ ] Avisar os autores ou registrar no texto: a matriz de MPD não é sempre euclidiana,
      ao contrário do que o artigo afirma (exemplo em `?phylo_beta`), e o picante não traz
      as funções que o artigo anuncia.

## v0.3.0: fluxo de estudos de simulação

Objetivo: reduzir o código que o usuário repete. As vignettes precisaram de laços e
funções escritos à mão, o que indica o que falta no pacote.

- [ ] Função para rodar grades de cenários com réplicas pareadas por árvore e
      paisagem e semente controlada, com paralelismo pelo Futureverse (`future`,
      `furrr` ou `future.apply` em `Suggests`).
- [ ] Diagnóstico de equilíbrio. Com competição filogenética e dispersão alta, a
      exclusão regional continuou até 3000 passos, e hoje nada avisa o usuário.
- [ ] Saída em formato longo (tibble) a partir de `history`.
- [ ] Benchmark de `mc_dynamics()`. Avaliar código compilado só se o benchmark
      mostrar ganho real, como foi feito no treeshapesim. O ABC depende disso.
- [ ] `summary_stats()`: as mesmas estatísticas para dados simulados e observados
      (alfa, beta, gama, ocupação, SAD, equitabilidade, SES MPD e MNTD, decaimento com
      a distância, relação entre traits e ambiente) e, com mais de um tempo, variação
      temporal, como o CV das abundâncias. Em Guzman et al. (2022), estatísticas
      temporais aumentaram o poder explicativo em até 59%, e a dispersão dependeu delas.
- [ ] Piloto de identificabilidade, com poucos parâmetros: quanto dá para recuperar a
      partir dessas estatísticas, com um levantamento único e com séries temporais? É o
      ponto de decisão para o ABC, e o desenho deve ser combinado antes de rodar.

**Concluído quando:** as vignettes puderem ser reescritas com essas funções, com bem
menos código.

## v0.4.0: base do hub

Objetivo: uma interface comum, para que simulação e métodos troquem dados entre si e
com outros pacotes. As versões seguintes dependem desta.

- [ ] Estrutura de dados comum: comunidade, árvore, traits, ambiente e coordenadas,
      com checagem de nomes. Serve para dados simulados e reais.
- [ ] Árvore e distribuições como entrada: pool regional por sítio e estado inicial
      vindos de fora. Cobre gen3sis, modelos de SSE e reconstruções biogeográficas.
- [ ] `species_traits()` aceitar traits vindos de fora, por exemplo do mvMORPH.
- [ ] Registrar a história dos traits ao longo dos ramos, se hoje só houver as pontas.
- [ ] Integração com o Herodotools (Nakamura et al. 2023): saídas nos formatos que ele
      lê e as métricas de `calc_tip_based_trait_evo()` como resumo, com vignette.
      Com a história conhecida, a simulação pode validar esses métodos.

**Concluído quando:** um cenário simulado passar pelo Herodotools e recuperar a
história conhecida.

## v0.5.0: métodos de metacommunity phylogenetics

Objetivo: levar para o R os métodos que só existem em MATLAB. O código está nos
suplementos de Peres-Neto et al. (2012) e Leibold et al. (2010), e nenhum dos dois
métodos está em pacote de R.

Vem antes do ABC por decisão de 2026-11-06. O código está em
`~/Downloads/ECOL_93_175/` (Peres-Neto et al. 2012) e em
`~/Downloads/ele_1523_sm_datas2.doc` (Leibold et al. 2010), fora do repositório porque
os suplementos não trazem licença. As comparações usaram Octave 11.3, instalado nesta
máquina.

- [x] Peres-Neto et al. (2012): `phylo_vector()` e `trait_partition()`, conferidas no Octave e nas Tabelas 1
      a 4 do artigo.
- [x] Tabela 5 com `logistic_gradients()` (Apêndice C), em `data-raw/validation_table5.R`.
      Frações com os vetores que geraram os dados e rejeição de E e S pelo teste global de
      `phylo_vector()`: o padrão entre cenários se repete, e I, III e IV ficam mais perto. A
      razão de variância não bate (0,62 contra 0,39 a 0,94), e o cenário VII tem `b` maior.
      A diferença está no artigo de validação; para fechar, só perguntando aos autores.
- [x] Semiparciais de `trait_partition()`, do Apêndice A, com testes. O suplemento não
      tem código para elas, então a conferência foi feita contra regressões separadas.
- [ ] `node_test()` roda nós em paralelo com `cores`. Se ainda for lento em árvores
      grandes, vetorizar `node_fit()`.
- [x] Leibold et al. (2010): `node_test()` e `node_contrasts()`. R2, R2 ajustado e F iguais
      ao MATLAB no Octave 11.3 (diferença abaixo de 1e-15), e os 13 contrastes de
      `CreateNodeContrasts` iguais. Em `data-raw/validation_node.R` (refeito na v0.5.2, com F), por nó, sites .053 e
      espécies .049 sob a nula dupla; p máximo .049 só com ambiente e .021 só com filogenia.
      Poder por nó .07, e ao menos um nó em 76% das réplicas, mas 24% sob a nula: é preciso
      corrigir para testes múltiplos. Os erros do original foram confirmados na leitura e
      corrigidos, com a diferença descrita em `?node_test`.
- [x] Semiparcial, FDR e eventos biogeográficos conforme o texto de Leibold et al. (2010).
      O MATLAB multiplica pelo módulo da inclinação e usa preditores brutos; o artigo
      divide e usa padronizados, e seguimos o artigo.
- [x] Codificação com comprimentos de ramo, na seção III do suplemento de Leibold et al.:
      `branch_lengths = TRUE`, conferida na Figura S2.
- [x] `plot_node_test()`: Figura 3 de Leibold et al. (2010), em filograma ou, para árvores
      muito grandes, leque, com os nomes das variáveis ao longo dos raios, com
      `ape::plot.phylo()`, sem dependência nova. O artigo não define a família do FDR:
      corrigimos preditores e eventos separadamente, o que está no help.
- [x] `plot_trait_partition()`: frações em barras, com negativas como zero (Legendre 2008).
      A tabela mantém os valores negativos.
- [x] `node_test()`: refúgios, 9999 permutações por padrão, permutação conjunta do artigo
      (`p_joint`, que rejeita 76% dos nós só com ambiente, contra 4,8% de `p_max`).
- [x] `daphniids`: dados de Leibold et al. (2010), refeitos da pesquisa pública da EPA.
- [x] Conferência de espécies por nome nas funções com árvore, traits ou contrastes.
- [x] Saíram `fourth_corner_test()` (fica no `ade4`) e `spatial_eigenvectors()` (fica no
      `adespatial`). `trait_partition()` cita ter Braak (2019) e o TraitEnvMLMWA para o
      teste de cada traço.
- [x] `plot_node_heatmap()`: nós por preditores com `semipartial`, como a Tabela S1.
- [ ] Estilo tidyverse na entrada e na saída, com álgebra de matrizes por dentro;
      permutações em paralelo com o Futureverse.
- [x] Testar a tradução contra o MATLAB em dados fixos, antes de corrigir qualquer
      coisa: R e Octave concordam até 1e-13, e as divergências estão nos helps.
- [x] Artigo "Analysing clades and traits in a metacommunity", com dafinídeos e simulação.

Situação em 2026-11-06: marco fechado na 0.5.1. Sobra o estilo tidyverse e o Futureverse,
se ainda fizerem sentido.

**Concluído quando:** a tradução reproduzir os resultados do MATLAB e os testes
reproduzirem o erro tipo I publicado.

## v0.5.2: biogeografia em `node_test()` e auditoria

Objetivo: tirar a maior barreira de uso do método de Leibold et al. (2010), que é montar a
matriz de eventos, e revisar o pacote do ponto de vista de quem usa. Vem antes de v0.3.0 e
v0.6.0 por decisão de 2026-11-06.

- [x] `biogeography_matrix()`: matriz de eventos a partir da área de cada site ou de uma
      matriz 0/1 de sites por áreas. Pares de áreas definidos a priori são o caminho
      principal, como em Leibold et al. (2010) e Leibold & Chase (2018, cap. 9), que usam
      hipóteses publicadas. "Cada área contra as outras" fica como opção exploratória, com
      fora das áreas como 0 (script de 2022) ou como -1 (script de 2015), e aviso quando um
      lado tiver poucos sites. Os dois scripts estão em `legacy/` do af-frogs-climate-rarity.
      Sobrepor polígonos e grade fica para o usuário, com `sf::st_join()`, sem dependência
      nova. Reproduz os seis eventos de `daphniids` e a tabela de refúgios de 2022 (com 5 como 1).
- [x] `site_permutations` em `node_test()`: para eventos, cada linha move as áreas pelos
      sites. Em `data-raw/validation_events.R` (1000 réplicas, toro de 15 por 15, blocos de 25
      sites sem efeito, comunidade seguindo ambiente oculto autocorrelacionado), rejeição por
      nó a 5%: .451 com permutação livre e .048 com o toro; ao menos um nó em 89% e 26% das
      réplicas, sem correção. Só testado com toro.
- [x] Nulo restrito para paisagens com borda: `patch_permutations()` e `area_permutations`
      em `node_test()`. Nenhum método publicado cobre permutação exata, bordas, área
      irregular e variável binária, então é método próprio. Em
      `data-raw/validation_edges.R` (200 réplicas; quadrado, faixa e banda irregular; blocos,
      áreas em manchas e área contra as outras), rejeição por nó a 5% sem efeito: .030 a .063,
      contra .44 a .55 com permutação livre. Descartados no protótipo: curva de Hilbert com
      deslocamento cíclico (.16 a .19), rotação rígida (lenta; .06 a .12) e mancha única por
      lado (conservador em áreas com várias manchas, .02 a .03, e menos poder). Falta aplicar a
      dados reais.
- [x] `plot_node_test()`: triângulos para eventos, como a Fig. 9.5 de Leibold & Chase
      (2018), além da cor.
- [x] Citar Leibold & Chase (2018) em `biogeography_matrix()`.
- [x] Auditoria do código com foco em experiência de uso, com dois revisores (2026-11-06).
      Corrigidos os erros que mudam resultados e as falhas silenciosas: F nos sorteios de
      `node_test()` (a estatística parcial do MATLAB dependia da escala), nós com clado
      ausente, espécies ausentes em `pcc_test()`, `phylo_beta_test()` e `trait_partition()`,
      `cores`, sites por nome, checagens em `R/checks.R` e `permutations` nos resultados.
- [x] `plot_trait_partition()`: barras sobrepostas, como a Fig. 1 de Desdevises et al. (2003).
- [x] Exemplos, README e vinhetas com `treeshapesim::sim_tree_shape()`. Testes seguem com
      `ape`. Os `.rds` das vinhetas foram refeitos: com a árvore nova, SES MPD fica em -1,3
      (era -1,7), e o toro rejeitou .08 e .07 (eram .035 e .050, com erro padrão .015).
- [x] Nota "Adapted from the MATLAB" em `node_test()`, `node_contrasts()`, `phylo_vector()` e
      `trait_partition()`. Ver a pendência sobre licença do código.
- [x] Barra de progresso com `purrr` e estimativa de tempo em sessões interativas, em
      `node_test()`, `pcc_test()`, `phylo_beta_test()`, `phylo_vector()` e `trait_partition()`.
      `purrr` entra em `Imports`.

Da auditoria, para antes da 1.0:

- [x] Descoberta: tabela de pergunta para função no README, página `?comdynamite`, novo Title
      e Description, artigo sobre preparar dados próprios e artigo de `pcc_test()` e
      `phylo_beta_test()`.
- [x] Nomes: `level` no lugar de `alpha` para significância, `n_perm` com padrão 999, `env` em
      `phylo_vector()`, `p_species`, `p_partial_max`, `r2_mean` em `pcc_test()` e `r2_adjusted`
      em `phylo_beta_test()`. `statistic` fica: é a estatística de cada teste, descrita na ajuda.
- [x] Métodos `summary()`, `plot()` e `as.data.frame()` para os resultados, em `R/methods.R`.
- [x] R2 de `node_test()` acima de 1 (herdado do original): aviso na ajuda e `r2_weighted`.
- [x] `site_permutations` em `trait_partition()`, e sites vazios com permutações restritas em
      `pcc_test()`, `phylo_beta_test()`, `phylo_vector()` e `trait_partition()`.
- [x] Baixa gravidade: ramos zero em `node_contrasts()`, nomes repetidos e heatmap só com NA,
      `logistic_gradients()` com poucos sites, `mc_dynamics()` com um site, `r2` NaN em
      `pcc_test()`, `n_draws`, erro dentro de `cores`, checagens do simulador.
- [x] Formas dos traits, na auditoria com Fable (2026-11-07): `alpha` alinhado por nome em
      `simulate_metacommunity()`, `as_traits()` para matriz e tibble, aviso de nomes inventados.

Situação em 2026-11-06: marco fechado na 0.5.2. O nulo com borda e a descoberta saíram na
0.6.0.

Opção, sem data: calanoides como segundo exemplo, se as espécies estiverem nas tabelas
públicas da EMAP. A topologia sairia da Fig. 9.5b, e os resultados do livro não foram
publicados em artigo, então seria exemplo novo, sem afirmar que reproduz o livro. O texto
do livro diz que os dafinídeos não mostram restrição biogeográfica, mas a Fig. 9.5a tem
vários triângulos: vale perguntar aos autores.

**Concluído quando:** a matriz dos refúgios e a dos dafinídeos saírem da função, o erro tipo
I dos eventos estiver medido com sites agregados no espaço e os achados da auditoria
estiverem resolvidos ou registrados aqui.

## v0.6.0: nulo com borda e descoberta (concluída em 2026-11-06)

`patch_permutations()`, `area_permutations` em `node_test()`, novo título e documentação por
pergunta. Os detalhes estão na seção da v0.5.2, onde os itens nasceram.

- [x] `patch_permutations()` com autocorrelação na escala da paisagem (bandas de 6, 12 e 25
      sites, 200 réplicas, sem efeito): rejeição de .022 a .060, contra .49 a .75 com
      permutação livre e .056 a .095 com o toro no quadrado. Conservador para área contra as
      outras (.022 a .039). Em `data-raw/validation_edges_scale.rds`.
- [x] Aplicado aos refúgios da Mata Atlântica, fora deste repositório (dados não publicados).
      Com células esparsas, o resultado depende muito de `gap`.
- [x] `patch_permutations()` não preserva a distância entre os lados do evento: resolvido em
      parte com `rigid_permutations()` (movimento rígido do evento inteiro; `independent` para
      lados distantes). Lados próximos: .022 a .048. Eventos que ocupam a extensão da paisagem
      ficam sem layouts ou, com `independent`, liberais (.06 a .09; .12 com autocorrelação longa).
      Não há nulo válido para lados distantes sem conhecer o alcance da autocorrelação.
      Em `data-raw/validation_rigid.R`. Registro do problema: com áreas fixas
      no quadrado 20 por 20: adjacentes .004 a .02; cantos opostos e metades .06 a .10 (.065 e
      .064 na minha conferência, 150 réplicas; .078 a .101 na auditoria, 60 réplicas). A
      validação sorteava as áreas, o que tira a média disso. Candidato: mover o conjunto de
      manchas de um evento junto, com deslocamento e rotação comuns.
- [x] Ramo fora de grade de `rigid_permutations()`: agora enumera rotações (passo que move o
      sítio mais distante em cerca de `gap`, de 4 a 360), reflexões e translações para cada
      sítio, sobre um raster com fechamento (buracos contam como dentro, a borda não), em vez
      de sortear movimentos; o conjunto de encaixes é fixo. Em `data-raw/validation_offgrid.R`
      (150 réplicas, 99 permutações): áreas adjacentes .047 e .059; dois blocos aleatórios em
      quatro paisagens .068 a .080, liberal; `patch_permutations()` .012 a .063; poder 23%
      contra 12%. Fica pendente: blocos aleatórios fora de grade acima de .07, e eventos
      grandes em grade com jitter são lentos (225 sítios em 750: 701 encaixes em 12 min).
      Nos refúgios, a parte grande tem um único encaixe, de forma determinística. Conferência
      independente (região em L, `data-raw/validation_nulls_check.R`): .070 e .087 fora de
      grade, contra .055 e .066 com `patch_permutations()`; grade esparsa .057 e .061.
- [x] Critério para `gap`: mediana da distância ao quarto sítio mais próximo (passo da grade
      cheia; espaçamento com sítios esparsos). `patch_permutations()` usa 1,5 vez isso (cantos
      com 25% a 40% da grade: .056 a .059, contra .15 a .17 com o padrão antigo; holdout .021
      a .059). `rigid_permutations()` usa o valor direto e não erra mais com poucos buracos
      (.055 a .061 com 40% a 80% dos sítios). Em `data-raw/validation_gap.R`. Cantos opostos
      esparsos e clusters seguem sem encaixes para o rígido.
- [x] Testes das frações `a` e `c` de `trait_partition()`. O MATLAB permuta ambiente e espaço
      com a mesma ordem de sítios e compara a própria fração (calcula o F parcial e não usa).
      Em vez de trocar a permutação, trocamos a estatística pela parcela do componente que o
      outro conjunto deixa (monótona no F parcial), o que mantém `site_permutations` e a
      permutação de espécies. No tuning (200 réplicas, grade 8 por 8) ela equivale ao
      Freedman-Lane e à permutação só de `env`; resultados do holdout (500 réplicas) em
      `data-raw/validation_fractions.R` e no artigo de validação. Limites registrados: com
      `env` estruturado em escalas fora de `space` o teste livre de `a` infla, e o toro
      segura; o toro é quase degenerado para `c`, porque mover os MEMs rigidamente quase não
      muda o espectro deles.
- [x] Teste da fração compartilhada `b` com `shared_null`: `"msr"` (Moran spectral
      randomization de `env` a partir de `coords`, via `adespatial::msr()` com a base de
      `dbmem()`) ou `"sites"` (`env` sozinho movido por `site_permutations`). Só testa se `env`
      estiver relacionado a `space` (`level`), como `envspace.test()`. p de sítios do MSR, p de
      espécies da permutação de traits já existente, `p_max` o maior; a permutação de espécies
      responde se o trait importa, não se a estrutura compartilhada é espúria, e o teste de
      `p_max` precisa dos dois (com trait oculto, o MSR sozinho rejeita 54%). Holdout com
      `env` estruturado e sem efeito: permutação livre de `env` .48, MSR .056, toro .11
      (liberal, sem explicação; fica como limite). Poder: 1.0 só ambiente, .43 com espaço
      agindo também, .28 com ambiente fraco. Referência verificada: Bauman et al. (2019,
      Oikos); Clappe et al. (2018) trata da partição com MSR, não do teste de `b`.
- [ ] Nulo restrito para `c` em `trait_partition()`: `site_permutations` move `env` e `space`
      juntos, e mover MEMs de uma grade quase não muda o que explicam; no holdout o toro foi
      liberal para `c` com `env` forte. Candidato: mover só `env` nos testes de `a` e manter
      `c` com permutação livre, ou documentar `c` como descritivo com nulo restrito.

## Inferência por ABC (planejada como v0.6.0, número a definir)

Objetivo: estimar parâmetros e comparar modelos a partir de dados observados. Depende do
piloto da v0.3.0 e da estrutura comum da v0.4.0.

Ler antes de implementar: Guzman et al. (2022), lido em parte; mcfly (Duarte et al. 2023); CAMI (Ruffley et al.
2019), pacote de R com random forests e ABC para montagem de comunidades; van der Plas
et al. (2015). Definir o que o comdynamite acrescenta a esses trabalhos.

Pontarp et al. (2019), lido com o Appendix S4: propõem um ciclo
de construir o modelo, analisá-lo sem dados, estimar parâmetros e avaliar. O caso de
estudo usa ABC por rejeição com 14 mil simulações por modelo, distância euclidiana entre
estatísticas (riqueza, abundância, distribuição de traits, NRI, NTI) escaladas pelo desvio
padrão, priors amplos e dois modelos, escolhidos pela razão de aceitação, tudo em dados
simulados e só com a comunidade final. O limiar sai do gráfico das distâncias e de uma
análise de sensibilidade; vale copiar esse passo, mas ele muda entre os trechos (1 na
Fig. 5, 3 e 5 no apêndice). O que levar para cá:

- Se um processo não é recuperado em dados simulados, não deve ser inferido em dados
  reais. É o papel do piloto da v0.3.0 e da recuperação abaixo.
- Traits e filogenia juntos separaram processos com sinais parecidos, então
  `summary_stats()` deve cobrir os dois.
- O modelo errado deu estimativas aparentemente boas, com a largura de nicho
  superestimada. Uma boa estimativa não valida o modelo: daí os PPCs.
- Eles não selecionaram estatísticas e citam Prangle et al. (2014) para isso. A razão de
  aceitação depende do limiar e das estatísticas (Robert et al. 2011), então preferir
  `abcrf` e medir o erro de classificação em dados simulados.
- O que o comdynamite acrescenta: espaço com muitos sítios, pacote em R (o código deles
  é MATLAB), cobertura dos posteriores e checagem do modelo, que o caso de estudo não faz.

- [ ] Tabela de referência: sortear parâmetros dos priors e simular, com a árvore, o
      ambiente e as coordenadas observados e com detecção imperfeita.
- [ ] Estimação e escolha de modelos com `abc` (Csilléry et al. 2012) e `abcrf` (Pudlo
      et al. 2015), em `Suggests`. Random forests lidam bem com muitas estatísticas
      correlacionadas e ordenam a importância delas.
- [ ] Custo: cerca de 0,2 s por simulação pequena (36 espécies, 20 sítios, 300 passos),
      ou 6 h por 100 mil simulações em um núcleo. Guzman et al. (2022) usaram 100
      manchas, 50 espécies e 2200 passos.
- [ ] Recuperação em dados simulados com parâmetros conhecidos: viés e cobertura.

Avaliação do modelo, sempre em conjunto com a inferência:

- [ ] Posterior predictive checks, com estatísticas que não entraram no ajuste.
- [ ] LOO sobre a tabela de referência: `cv4abc()` e `cv4postpr()` do `abc`, com cada
      simulação como pseudo-observação; `gfit()` para adequação global.
- [ ] K-fold por sítios ou blocos espaciais: ajustar com parte dos sítios e prever as
      estatísticas dos demais. Blocos evitam que a autocorrelação espacial infle o acerto.
- [ ] Diagnóstico de Pareto k dos pesos de importância, com `loo::psis()`, para saber
      quando a posterior ponderada depende de poucas simulações.
- [ ] PSIS-LOO padrão (Vehtari et al. 2017) precisa de log-verossimilhança por observação,
      que o ABC não dá. Avaliar likelihood sintética (Wood 2010) por sítio ou por
      estatística, que daria esse termo de forma aproximada.

Limites a deixar claros: um bom ajuste não prova o mecanismo, e dados de um único
levantamento podem não informar a dispersão.

**Concluído quando:** uma vignette levar dados reais do ajuste à avaliação, com
cobertura medida em dados simulados.

## v0.7.0: evolução de traits e modelos nulos

Objetivo: ampliar os modelos de traits e os nulos com dinâmica histórica, chamando os
pacotes de referência em vez de reimplementá-los.

- [ ] Mais opções de evolução de traits, multivariadas e com covariância entre traits,
      pelo mvMORPH (Clavel et al. 2015).
- [ ] Funções que usam traits, como `assemble_community()`, com mais de um trait e com
      restrições na estrutura entre traits, como em `mvSIM()` e nos ajustes do mvMORPH:
      sem covariância (`"diagonal"`), taxas iguais (`"equal"`), variâncias ou correlações
      iguais (`"variance"`, `"correlation"`), matrizes proporcionais (`"proportional"`) ou
      compartilhadas (`"shared"`), e matriz definida pelo usuário. Definir como cada processo
      usa vários traits: distância multivariada na similaridade limitante, nicho
      multidimensional no filtro.
- [ ] Modelo nulo dinâmico do DAMOCLES (Pigot & Etienne 2015): colonização e extinção
      local ao longo da árvore.

**Concluído quando:** cada opção tiver teste e exemplo em vignette.

## v0.8.0: traits, filogenia e biogeografia

Objetivo: transformar o desenho de Benício et al. (2021) em um fluxo do pacote e testar,
com a verdade conhecida, se ele separa tempo para especiação, taxas de diversificação
e conservação de nicho. Scripts e dados do artigo estão no Dryad
(https://doi.org/10.5061/dryad.sf7m0cg65).

Simular o cenário:

- [ ] Árvore e regiões geradas juntas, por exemplo com um simulador de GeoSSE, em cenários
      que variem especiação, extinção, dispersão e tempo de colonização por região.
- [ ] Traits discretos e filtragem por eles, como os modos reprodutivos. Hoje o modelo
      só tem traits contínuos.
- [ ] Duas regiões com pools próprios e dispersão limitada entre elas, usando a
      entrada da v0.4.0.
- [ ] Amostragem imperfeita: sítios, dias de coleta e espécies não detectadas.

Analisar, chamando os pacotes de referência sobre a estrutura comum:

- [ ] Riqueza por rarefação e extrapolação (`iNEXT`).
- [ ] Sinal filogenético em traits discretos: Rao decomposto nos nós (`adiv`, `ade4`).
- [ ] Diversificação por região: GeoHiSSE (Caetano et al. 2018), no lugar do GeoSSE
      usado no artigo, e testes baseados nas pontas, como FiSSE (Rabosky & Goldberg
      2017) e ES-sim (Harvey & Rabosky 2017).
- [ ] Curvas de linhagens no tempo por região (`ape`).
- [ ] DEC (BioGeoBEARS; Matzke 2014) e tempo de colonização por idades crown e stem.
- [ ] MPD, MNTD e análogos funcionais, com o pool definido em três extensões: só uma
      região, só a outra e as duas juntas (`picante`).
- [ ] Relação da diversidade com clima, clima passado e topografia nas três extensões.
- [ ] Diversificação em função de clima e topografia, com o RPANDA (Morlon et al. 2016).

Depois:

- [ ] Avaliar o fluxo: em que cenários ele recupera a causa simulada e com que frequência
      erra. Incluir árvore incompleta e poucos sítios, limitações citadas no artigo, e
      os falsos positivos de modelos SSE (Rabosky & Goldberg 2015). Comparar GeoSSE,
      GeoHiSSE e testes baseados nas pontas nos mesmos cenários.
- [ ] ClaSSE (Goldberg & Igić 2012) e outros modelos de Markov com k estados.
- [ ] SBEARS, do Herodotools. Conferir o que ele faz: entrou aqui como modelo contínuo,
      e o nome é Site Based Estimation of Ancestral Range of Species.

**Concluído quando:** um estudo de simulação mostrar quando o fluxo acerta e quando erra.
Esse estudo pode virar artigo próprio.

## v0.9.0: mais processos

Objetivo: ampliar os cenários do simulador. Escolher pelo uso real.

- [ ] Vários eixos ambientais (hoje há um só).
- [ ] Paisagens irregulares ou lidas de coordenadas, além da grade.
- [ ] Null models restritos com bordas para preditores contínuos e para o simulador.
      `patch_permutations()` cobre só eventos; candidato: Constant Moran permutation
      (Arthur 2024), a testar.
- [ ] Gerador de ambiente temporal fiel ao de Thompson et al. (2020), como alternativa
      a `env_over_time()`.

**Concluído quando:** os itens escolhidos tiverem testes e exemplo em vignette.

## v1.0.0: CRAN e artigo

Objetivo: API estável e pacote no CRAN.

- [ ] Congelar a API e registrar as funções depreciadas no `NEWS.md`.
- [ ] Tirar `Remotes` do `DESCRIPTION`: o CRAN não aceita. `treeshapesim` precisa estar
      no CRAN, ou ficar só em `Additional_repositories` com o CI ajustado.
- [ ] `cran-comments.md`; rodar `R CMD check --as-cran`, `rhub` e `win-builder`.
- [ ] Conferir o tempo dos exemplos e das vignettes.
- [ ] Submeter ao CRAN e criar a tag `v1.0.0`.
- [ ] Artigo de software. Periódico e momento a definir; depende do manuscrito do
      pcps-sim (ver abaixo).

**Concluído quando:** o pacote estiver aceito no CRAN.

## Pendências e decisões

- **Imigração na grade.** Decidir se `mc_dynamics()` passa a conservar dispersores, como o
  mcomsimr, ou mantém a Poisson do código original. Mudar altera os resultados de
  sementes antigas e exige refazer as vignettes.
- **Erros no artigo de Thompson et al. (2020).** Tabela 1 (beta), faixa de dispersão,
  conversão do ambiente em postos e metade dos imigrantes no cenário de
  competição-colonização não batem com o código. Vale citar isso no artigo do pacote.
- **Dependências do hub.** Hoje o pacote só importa `ape`, `purrr` e pacotes base (`graphics`, `grDevices`, `parallel`, `stats`). Tidyverse, Futureverse
  e pacotes de biogeografia pesam em `Imports`; preferir `Suggests` sempre que der.
- **Código em MATLAB.** Os suplementos não trazem licença. Reimplementar a partir das
  equações e citar; se for copiar a estrutura do código, falar com os autores.
- **Protocolo estático.** O texto de Peres-Neto et al. (2012) dá `h` entre 0,5 e 1 e usa
  **L** como abundância, sem Poisson nos trechos lidos (a etapa 1 não veio). O pcps-sim
  usa `h` entre 0 e 30 e Poisson, e o manuscrito atribui isso a Peres-Neto et al.
  (2012): conferir essa citação, talvez a fonte seja Dray & Legendre (2008).
- **Appendix C.** Não dá a distribuição de `a` nem a escala de `b`; diz que a presença
  vem de um valor "normally generated", o que gera muito ruído (a legenda da Fig. 4
  comenta o ruído). `simulate_logistic()` usa uniforme por padrão e aceita
  `draw = "normal"`. O cenário III é ambíguo. O código de Ives & Helmus (2011,
  `PGLMM_model_simulations.m`, em `~/Downloads/3567789/`) esclarece parte disso: `a` é
  N(0; 0,5), `b` tem desvio 1 na amostra de espécies, a presença vem de `rand` (uniforme),
  e a redução U(0,5; 1) entra uma vez por gradiente, enquanto o Appendix C a aplica uma
  vez. `simulate_logistic()` adotou esses valores na 0.5.1, e a redução fica uma só, como
  no Appendix C. O resultado da Tabela 5 está no artigo de validação.
- **treeshapesim.** O artigo de integração com o comdynamite usa `species_traits()`,
  `competition_matrix()`, `landscape_implicit()` e `simulate_metacommunity()`: mudar
  assinatura ou saída delas quebra o site do treeshapesim. Quando `prune_tips()` estiver
  na `main` de lá, pôr o link do artigo no README ou em `evaluating-methods`. O
  comdynamite fica só em `Config/Needs/website` no treeshapesim, sem dependência circular.
- **Sessão anterior.** A sessão que começou o pacote não foi avisada por aqui de que o
  desenvolvimento mudou de lugar. Commits novos dela vão divergir do `main`.
- **pcps-sim.** Se a submissão for duplo-cega, avaliar se o `NEWS.md` pode manter "Code
  from the pcps-sim project". A vignette `evaluating-methods` compara permutações
  livres e do toro, tema próximo do piloto de validade sob H0 daquele projeto.
- **Resultados pré-calculados.** Os `.rds` em `vignettes/` vêm de
  `data-raw/vignette_ecology.R` e `data-raw/vignette_method.R`. Se o modelo mudar,
  é preciso rodar os scripts de novo e revisar os números citados no texto.
- **Páginas órfãs do site.** O deploy do pkgdown usa `clean: false`, então a página de uma
  função removida fica no `gh-pages` e precisa ser apagada à mão.
- **`data-raw/hex_dynamite.R`.** Sobrescreve o `logo.png` oficial se for rodado. Remover
  ou mudar o arquivo de saída.

## Referências citadas neste roadmap

- Desdevises, Y., Legendre, P., Azouzi, L., & Morand, S. (2003). *Evolution*, 57, 2647–2652. https://doi.org/10.1111/j.0014-3820.2003.tb01508.x
- Leibold, M. A., & Chase, J. M. (2018). *Metacommunity Ecology*. Monographs in Population Biology 59. Princeton University Press. https://doi.org/10.2307/j.ctt1wf4d24

Conferidas no Crossref em 2026-11-05.

- Benício, R. A., et al. (2021). *Journal of Biogeography*, 48, 2211–2225. https://doi.org/10.1111/jbi.14145
- ter Braak, C. J. F. (2019). *Methods in Ecology and Evolution*, 10, 1962–1971. https://doi.org/10.1111/2041-210X.13278
- Caetano, D. S., et al. (2018). *Evolution*, 72, 2308–2324. https://doi.org/10.1111/evo.13602
- Csilléry, K., et al. (2012). *Methods in Ecology and Evolution*, 3, 475–479. https://doi.org/10.1111/j.2041-210x.2011.00179.x
- Dray, S., & Legendre, P. (2008). *Ecology*, 89, 3400–3412. https://doi.org/10.1890/08-0349.1
- Clavel, J., et al. (2015). *Methods in Ecology and Evolution*, 6, 1311–1319. https://doi.org/10.1111/2041-210x.12420
- Duarte, L. D. S., et al. (2016). *Methods in Ecology and Evolution*, 7, 937–946. https://doi.org/10.1111/2041-210x.12547
- Duarte, L. D. S., et al. (2023). *Ecography*, 2023. https://doi.org/10.1111/ecog.06544
- Goldberg, E. E., et al. (2011). *Systematic Biology*, 60, 451–465. https://doi.org/10.1093/sysbio/syr046
- Goldberg, E. E., & Igić, B. (2012). *Evolution*, 66, 3701–3709. https://doi.org/10.1111/j.1558-5646.2012.01730.x
- Guzman, L. M., et al. (2022). *Ecology*, 103. https://doi.org/10.1002/ecy.3683
- Harvey, M. G., & Rabosky, D. L. (2017). *Methods in Ecology and Evolution*, 9, 984–993. https://doi.org/10.1111/2041-210x.12949
- Ives, A. R., & Helmus, M. R. (2011). *Ecological Monographs*, 81, 511–525. https://doi.org/10.1890/10-1264.1
- Legendre, P. (2008). *Journal of Plant Ecology*, 1, 3–8. https://doi.org/10.1093/jpe/rtm001
- Leibold, M. A., et al. (2010). *Ecology Letters*, 13, 1290–1299. https://doi.org/10.1111/j.1461-0248.2010.01523.x
- Matzke, N. J. (2014). *Systematic Biology*, 63, 951–970. https://doi.org/10.1093/sysbio/syu056
- Morlon, H., et al. (2016). *Methods in Ecology and Evolution*, 7, 589–597. https://doi.org/10.1111/2041-210x.12526
- Nakamura, G., et al. (2023). *Journal of Biogeography*, 51, 533–543. https://doi.org/10.1111/jbi.14610
- Pontarp, M., et al. (2019). *Methods in Ecology and Evolution*, 10, 450–460. https://doi.org/10.1111/2041-210x.13129
- Prangle, D., et al. (2014). *Statistical Applications in Genetics and Molecular Biology*, 13, 67–82. https://doi.org/10.1515/sagmb-2013-0012
- Pudlo, P., et al. (2015). *Bioinformatics*, 32, 859–866. https://doi.org/10.1093/bioinformatics/btv684
- Rabosky, D. L., & Goldberg, E. E. (2015). *Systematic Biology*, 64, 340–355. https://doi.org/10.1093/sysbio/syu131
- Rabosky, D. L., & Goldberg, E. E. (2017). *Evolution*, 71, 1432–1442. https://doi.org/10.1111/evo.13227
- McArdle, B. H., & Anderson, M. J. (2001). *Ecology*, 82, 290–297. https://doi.org/10.1890/0012-9658(2001)082[0290:FMMTCD]2.0.CO;2
- Peres-Neto, P. R., & Kembel, S. W. (2015). *Plant Ecology*, 216, 709–724. https://doi.org/10.1007/s11258-014-0405-0
- Peres-Neto, P. R., et al. (2012). *Ecology*, 93, S14–S30. https://doi.org/10.1890/11-0494.1
- Pigot, A. L., & Etienne, R. S. (2015). *Ecology Letters*, 18, 153–163. https://doi.org/10.1111/ele.12395
- Robert, C. P., et al. (2011). *PNAS*, 108, 15112–15117. https://doi.org/10.1073/pnas.1102900108
- Ruffley, M., et al. (2019). *Ecology and Evolution*, 9, 13218–13230. https://doi.org/10.1002/ece3.5773
- van der Plas, F., et al. (2015). *Ecology*, 96, 1502–1515. https://doi.org/10.1890/14-0454.1
- Vehtari, A., et al. (2017). *Statistics and Computing*, 27, 1413–1432. https://doi.org/10.1007/s11222-016-9696-4
- Wood, S. N. (2010). *Nature*, 466, 1102–1104. https://doi.org/10.1038/nature09319
