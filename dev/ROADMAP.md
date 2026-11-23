# comdynamite: estado e próximos passos

Atualizado em 2026-11-22. Versão atual: 0.11.0 (tag `v0.11.0`, release, DOI
`10.5281/zenodo.24446713`), após a 0.10.0 (tag `v0.10.0`, release, DOI `10.5281/zenodo.24423302`) e a 0.9.0 (DOI
`10.5281/zenodo.24341231`; DOI de conceito `10.5281/zenodo.23993533`). Em desenvolvimento:
0.11.0.9000.
Tags: `v0.1.0`, `v0.2.0`, `v0.5.0`, `v0.5.1`, `v0.5.2`, `v0.6.0` a `v0.6.3`, `v0.7.0`, `v0.8.0`,
`v0.9.0`, `v0.10.0` e `v0.11.0`.
O marco v0.3.0 fechou em 2026-11-07 e saiu na 0.6.2, sem tag `v0.3.0`.

## Próximos passos (2026-11-22)

0. **0.11.0 em 2026-11-22** (tag `v0.11.0`): fecha os padrões do rOpenSci. Todos os padrões
   aplicáveis do `srr` documentados (blocos A a D, PRs #150 a #153), com as tags também nos
   arquivos de teste, e o pkgcheck local sem falhas nos checks principais (PR #155). Restam só
   avisos de estilo do goodpractice. Release no GitHub com DOI `10.5281/zenodo.24446713`. A
   submissão continua adiada (ver 1b).
   Depois da tag, na 0.11.0.9000: o workflow do pkgcheck no GitHub voltou a funcionar (PR #159,
   script próprio `.github/pkgcheck.R` na imagem da action, que contorna o erro do resumo de
   pkgcheck-action #77; roda à mão, uns 40 min) e o goodpractice passou nos checks de tamanho
   (até 50 linhas) e de complexidade (até 15) das funções (PR #158). Cerca de 105 funções em
   43 arquivos foram divididas em auxiliares internas, sem mudança de comportamento (resultado,
   mensagens e sequência de sorteios idênticos em comparações golden), e os testes seguem o
   estilo do lintr (888 `expect_identical()`). Ficam de propósito: as 307 comparações de ponto
   flutuante em `expect_equal()`, as linhas acima de 80 colunas (o pacote usa 100) e
   `shift_tree_once()` cerca de 25% mais lenta (interna, só em `data-raw`).
1. **0.10.0 em 2026-11-21** (tag `v0.10.0`): fecha o bloco A (nulos com borda e
   `env_surrogates()`, gerador de Thompson et al. 2020, competição por traits, nicho com eixos
   correlacionados, dispersão no tempo, tempo desde a mudança, ambiente sazonal e redes
   fluviais, dispersão e tolerância a distúrbio por espécie ou grupo) e traz o artigo de
   persistência em redes fluviais com os cenários do Pantanal. Release no GitHub com DOI
   `10.5281/zenodo.24423302`. As duas pendências do bloco A entraram antes da tag: nicho
   correlacionado com `evolution` e com especiação (item 4) e `kernel_group` no motor evolutivo
   e com seleção de imigrantes ou vários traços (item 8). Restam decisões menores: `sigma = 0`
   com correlação e a fissão por distância pelo kernel do grupo. Próximo marco de código: a
   escolher entre os blocos B a J.
   A 0.9.0 saiu em 2026-11-18 (tag `v0.9.0`, release, DOI `10.5281/zenodo.24341231`).
1b. **Padrões do rOpenSci (auditoria de 2026-11-21).** O pacote entra pela via de software
   estatístico (categorias EA, RE e BS, 206 padrões a anotar com `srr`). Fase 1 (correções
   rápidas: CONTRIBUTING, código de conduta, ORCID, badge, WORDLIST, `.lintr`, pkgcheck no CI,
   falhas do pkgcheck) na 0.10.0.9000. Fase 2 (atribuição do código derivado de MATLAB) no PR #142, com o
   e-mail aos autores enviado em 2026-11-21; mensagens de citação dos artigos originais no
   PR #147, a pedido de um dos autores. Fase 3 (validação de entradas) no PR #143 e fase 4
   (suíte estendida de testes) no PR #144. Fase 5, anotação `srr`, começou em 2026-11-22
   (branch `srr-anotacao`). A consulta prévia aos editores fica para depois da decisão do JBI
   e de alguns meses de histórico estável: em 2026-11-22 o rOpenSci pôs o lagprotest em
   espera (ropensci/software-review#829) por ter o repositório novo demais, e o comdynamite
   tem o mesmo problema (primeiro commit em 2026-11-05).
   Fase 4 (testes estendidos, `tests/testthat/test-extended-*.R`, PR #144) achou quatro
   problemas, tratados no PR seguinte: (a) intervalos do ABC-RF aproximados (propriedade do
   método; ajuda corrigida, `min_node_size` novo); o critério de "recuperável" dos scripts
   `data-raw/abc_*.R` e dos pilotos perdeu o teto de 95% de cobertura (R2 >= 0,5 e cobertura
   >= 85%, coluna `conservative` acima de 95%), aplicado aos rds salvos sem rodar de novo por
   `data-raw/abc_criterion_update.R`: 26 linhas passaram a recuperáveis, todas com R2 >= 0,64;
   o artigo do JBI descreve o critério antigo e já diz que ele falhou por excesso de
   cobertura, então só mudar o texto se a revisão pedir; (b) teste da média de
   `tstats_summary()`: corrigido (nulos de T_IC.IR e T_PC.PR permutam a região inteira); (c)
   T_PC.PR: o nulo publicado dá ~80% de sítios "altos" falsos quando o tamanho das populações
   varia entre sítios; `null_pc = "populations"` evita isso e virou o padrão (decidido em 2026-11-21);
   `community_index()` sorteava cada sítio da região de forma independente: os nulos regionais
   agora permutam a região inteira (a distribuição de cada sítio não muda); (d)
   `trait_flex()`: propriedade do método (atenuação e teste conservador com turnover), ajuda
   corrigida. Também resolvidos: sinal dos autovetores de `phylo_vector()` e
   `phylo_eigenvectors()` fixado (maior valor absoluto positivo; o sinal trocava no Windows) e
   exemplo de `damocles_partition_test()` reduzido (o custo está no ajuste; acelerar
   `damocles_prune()` fica como ideia).
2. **Preprint reescrito para a 0.9.0** (PR #112, mesclado em 2026-11-18, cb65a09): preprint
   postado e artigo submetido ao Journal of Biogeography (Methods and Tools, double blind) em
   2026-11-18. Os arquivos renderizados ficam no Box (`Meus Artigos/2026 - comdynamite JBI`),
   fora do git. Falta registrar o DOI do preprint em `CITATION` e no README quando sair.
   As limitações resolvíveis do manuscrito entraram nos blocos A, B, C e H (ver "Limitações do
   manuscrito" abaixo), para que a revisão já encontre parte delas feitas.
3. **DAMOCLES:** terminar a validação S3 com o ClaDS2 e escrever o plano do artigo
   (`dev/plano-artigo-damocles.md`). O proxy de persistência por traits (ideia do Januário)
   foi testado com TetrapodTraits, longevidade e ocupação dinâmica no NAAMP (PR #120,
   2026-11-19); nota enviada ao Januário em 2026-11-19, aguardando a reunião.
4. **Delatorre et al. (2015):** experimento virtual feito (local, `dev/delatorre2015/`, fora
   do git). Sorting no ambiente medido é descartado; sorting em eixo não medido contra
   neutralidade fica em aberto com as estatísticas do artigo, e a matriz de espécies
   resolveria (correlação entre espécies e primeiro eixo dos resíduos). Matriz pedida aos
   autores em 2026-11-17.
5. **Fase 4 da especiação:** comparação com o gen3sis e o MESS, só em `data-raw`.
6. **Pendentes de decisão do usuário:** issue #56 do RPANDA e e-mail à Melissa Guzman (não
   enviado); nome final e modos de `simulate_diversification()`; sugestões do ECY25-0086.

Fechado sem continuação: o simulador sazonal com passo mensal (issue #97). O piloto de
identificabilidade (`data-raw/seasonal_pilot/`) mostrou que dados mensais de uma paisagem de
13 poças mal separam os modelos de comunidade; issue fechada como "not planned" em 2026-11-17.
As versões abaixo são uma sugestão: cada marco tem um critério de conclusão, e a
ordem pode mudar. Tag no formato `vX.Y.Z`, com release no GitHub.

## Itens abertos e plano de implementação (levantamento de 2026-11-17)

Levantamento de todos os itens abertos do roadmap: as caixas não marcadas e os "Falta:" que
ficaram dentro de itens marcados como feitos, sobretudo na v0.7.0. Conferido contra o código
na `main` (f3aab10). As fases 1 a 3 da especiação e a G multivariada estavam com a caixa
aberta, mas foram feitas (0.8.0); foram marcadas. Os rótulos dos blocos não são números de
versão: o número sai quando o bloco fecha.

**Bloco A, nulos e ambiente (fechado na 0.10.0 em 2026-11-21).** Itens pequenos,
dentro do pacote, que melhoram o núcleo antes de congelar a API.

1. *Nulos restritos com bordas para preditores contínuos e para o simulador.* Hoje o toro numa
   grade com bordas rejeita 6% a 25% sem relação real em `pcps_test()`, e
   `patch_permutations()` cobre só eventos (áreas). Plano: (i) estudo em
   `data-raw/validation_continuous_nulls.R` com preditores contínuos autocorrelados em grade
   com bordas, paisagem irregular (`landscape_points()`) e região em L, comparando toro,
   `rigid_permutations()`, redes de Moran (MSR, `adespatial::msr()`, em `Suggests`) e a
   Constant Moran permutation (Arthur 2024; referência não verificada, conferir antes de
   citar); (ii) implementar só o que ficar dentro de 0,04 a 0,06 de rejeição, como opção
   `null = ` em `pcps_test()`, `node_test()` e `phylo_beta_test()`; (iii) para o simulador,
   uma função que gere ambientes nulos com a mesma autocorrelação (substitutos espectrais ou
   MSR) para testar se um padrão simulado depende do ambiente. Pendências antigas que entram
   aqui: blocos aleatórios fora de grade liberais (acima de 0,07) e eventos grandes em grade
   com jitter lentos (v0.6.0), e aplicar os nulos com borda a dados reais (v0.5.2).
2. *Gerador de ambiente temporal fiel ao de Thompson et al. (2020).* `env_over_time()` é um
   AR(1) independente por sítio, sem estrutura espacial no tempo. Plano: (i) ler como o
   ambiente é gerado no código Julia dos autores e em `mcomsimr::env_generate()` (estrutura
   espacial, autocorrelação temporal, conversão em postos, que a seção de erros do artigo já
   registra) antes de escrever qualquer coisa; (ii) `env_thompson()` (nome provisório) que
   devolva a matriz sítio × passo aceita por `env_time`, com semente reprodutível e sem
   `RandomFields` (fora do CRAN), por exemplo campo Gaussiano por decomposição espectral;
   (iii) validação contra o mcomsimr com o ambiente mudando no tempo, que é o "Falta um
   ambiente que mude no tempo" da comparação da v0.2.0, com os mesmos testes KS das 80
   comparações; (iv) usar o gerador na vignette de validação no lugar do ambiente constante.
   Feito em 2026-11-19 (PR #124): `env_space_time()`, validada contra o mcomsimr e usada na
   vignette de validação; GIF em `data-raw/animation_env_space_time.R`, no artigo de validação.
3. *Competição a partir da distância entre traits*: `competition_matrix(type = "trait")`, com
   o mesmo núcleo de sobreposição de nicho que o motor evolutivo já usa (v0.7.0, "Falta").
   Feito em 2026-11-19 (branch `bloco-a-itens-3-4`), com vários eixos pelo produto das
   sobreposições.
4. *Nicho com correlação entre eixos* (matriz de largura por espécie) em `growth_rates()`, e
   paisagens que gerem vários eixos ambientais de uma vez em `landscape_grid()` (v0.7.0).
   Feito em 2026-11-19 (mesma branch): `niche_correlation()` (uma matriz para todas as
   espécies ou uma por espécie) e `landscape_grid(n_axes, axis_cor)`, com `env_axes`. Feito em
   2026-11-21 (branch `versao-0.10.0`): nicho correlacionado com `evolution` e com
   especiação, com Omega = D R D no crescimento, na resposta à seleção, na seleção de
   imigrantes, na sobreposição de nicho (`niche_overlap` e `competition_matrix(type =
   "trait")`) e em `maladaptation()`; as filhas herdam a correlação da mãe, e R = identidade
   reproduz os eixos independentes bit a bit. Falta: `sigma = 0` com correlação (recusado) e
   `mutation_var` com vários traits (já recusado antes).
5. *Dispersão variável no tempo e filtragem de imigrantes com especiação e vários traits*:
   `dispersal_time`, `kernels` e `kernel_step` com `evolution`, e `immigrant_selection` com
   `speciation` e com G multivariada, hoje recusados (artigo "Dispersal and disturbance through
   time"). É a terceira limitação do manuscrito, a única que é só engenharia. Feito em
   2026-11-19 (PR #116, d6df7f1): `simulate_evolving_metacommunity()` aceita `dispersal_time`,
   `kernels` e `kernel_step`, e `immigrant_selection` funciona com especiação e com G
   multivariada; só `mutation_var` continua recusado no motor de especiação. Sem especiação,
   a dispersão variável no tempo com vários traits já funcionava, ao contrário do que o
   manuscrito diz; corrigir o texto na revisão.
6. *Saída de "tempo desde a mudança"* no motor evolutivo, a única sugestão do ECY25-0086 feita
   só em parte. Feito em 2026-11-20: `evolution$since_change` (fim e `record_steps`) e a
   coluna `since_change` de `maladaptation()`.
7. *Ambiente sazonal e correlacionado ao longo de redes fluviais* (pedido em 2026-11-19).
   Objetivo: investigar dinâmica e persistência de longo prazo em metacomunidades fluviais,
   como as de Altermatt (2013, Aquatic Ecology, doi:10.1007/s10452-013-9450-3) e Carrara et
   al. (2012, PNAS, doi:10.1073/pnas.1119651109). A dispersão em rede dendrítica já cabe em
   `landscape_matrix()`; falta o ambiente. `env_space_time()` é estacionário, sem ciclo, e usa
   distância euclidiana, de modo que afluentes próximos em linha reta mudam juntos. Plano:
   (i) [feito em 2026-11-20: `env_seasonal()`, `river_floods()` e `disturbance` como matriz sítios x passos] componente sazonal, um ciclo determinístico com amplitude por sítio (por exemplo,
   cota ou distância ao canal) somado a variação entre anos, compondo com
   `env_space_time(ranks = FALSE)` e `env_pulse()` para cheias; (ii) covariância ao longo da
   rede, pelos modelos tail-up de Ver Hoef et al. (2006, Environmental and Ecological
   Statistics, doi:10.1007/s10651-006-0022-8) e Ver Hoef & Peterson (2010, JASA,
   doi:10.1198/jasa.2009.ap08248), já que a exponencial com distância de rede não é sempre
   positiva definida (correção em 2026-11-20, lendo Peterson & Ver Hoef 2010, Ecology, e Som
   et al. 2014, Environmetrics, pelo conector da Wiley: a exponencial da distância ao longo do
   canal é o modelo tail-down exponencial, que é válido; o que falha com distância de rede são
   modelos como o esférico. Feito em 2026-11-20: `river_covariance()` com tail-up, tail-down,
   euclidiano e pepita como componentes de variância, e `env_river()`; os PDFs de Ver Hoef et
   al. 2006 (eq. 7) e Ver Hoef & Peterson 2010 (eq. 4 e seção 2.2.1) foram lidos depois, em
   2026-11-20, e confirmam o peso como raiz do produto das frações a montante); (iii) um gerador de redes dendríticas, que entra no item (decisão do usuário em
   2026-11-19), para variar o tamanho e a forma da rede com réplicas; candidatos: optimal
   channel networks e árvores aleatórias com ordem de Strahler controlada; ler antes como os
   trabalhos de Altermatt e Carrara geram as redes; deve devolver coordenadas, a matriz de
   distâncias ao longo do canal e a direção do fluxo, para alimentar `landscape_matrix()` e a
   covariância de (ii), aceitando também redes reais. (iii) feito em 2026-11-20: `river_network()` (crescimento aleatório com `branching`, ordem de Strahler e magnitude de Shreve) e `landscape_river()` (kernel com escalas a jusante e a montante, dispersão por terra opcional, aceita redes do OCNet ou reais; OCNet é GPL-3, fica fora das dependências). Os textos completos de Altermatt e Carrara não estavam acessíveis (PMC só com resumo); as referências entram só como contexto; (iv) estudo de persistência de longo prazo com
   cheias sazonais, com a dispersão a favor e contra a corrente. (iv) feito em 2026-11-21: artigo
   `vignettes/articles/river-persistence.Rmd`, com o estudo genérico (topologia, direção da
   dispersão, cheias) e cenários inspirados no Pantanal (pulso normal, seca de 12 anos,
   regulação; guildas com rotas e tolerância à secagem do item 8, calibradas por Pompeu &
   Godinho 2006; lagoas isoladas fora da cheia), com GIF e esquema do desenho. A lição do piloto sazonal
   (issue #97) vale aqui: o gerador serve para dinâmica e persistência, não para inferência
   com poucos sítios. Referências conferidas no OpenAlex em 2026-11-19; ler os artigos antes
   de citar resultados.
8. *Dispersão e tolerância a distúrbio por espécie ou grupo de espécies* (pedido em
   2026-11-20, a partir dos cenários do Pantanal do item 7d, em que as guildas de peixes só
   diferiam em quanto emigravam e a secagem de uma lagoa eliminava todas as populações). Feito
   em 2026-11-20: `kernel_group` com `kernel_step` como matriz passos x grupos (cada grupo de
   espécies segue o seu kernel em cada passo), `dispersal_site` (multiplicadores da emigração
   por sítio e passo) e `disturbance_sensitivity` (probabilidade de a população sumir quando o
   sítio é perturbado, uma por espécie). Os padrões reproduzem os resultados anteriores bit a
   bit. Feito em 2026-11-21: `kernel_group` com seleção de imigrantes e com vários traços
   (os imigrantes e os genes de cada espécie seguem o kernel do seu grupo) e no motor
   evolutivo (`simulate_evolving_metacommunity(kernel_group = )`: um grupo por espécie
   inicial, as filhas ficam no grupo da mãe, o fluxo gênico de coesão segue o kernel do grupo
   e a saída devolve o grupo de cada espécie). Com todos os grupos na mesma matriz, o
   resultado é idêntico bit a bit ao de um `kernel_step` vetorial. Falta: a fissão por
   distância continua usando o kernel de `landscape` para todos os grupos (como já fazia com
   `kernels`); decidir se deve seguir o kernel do grupo. Cenários do Pantanal refeitos com as
   guildas em 2026-11-21 (item 7, iv).

Ordem sugerida para o que falta: 2, 1, 3, 4, 6, 7. Os itens 2, 3 e 4 foram feitos em 2026-11-19; o estudo (i) do item 1 foi feito em 2026-11-20 (200 réplicas, seis nulos: livre, toro, MSR, shift-rotate ad hoc, substituto paramétrico por postos e variograma de Viladomat et al. 2014; a Constant Moran permutation saiu porque Arthur 2024/2025 trata de nulos de redes). Nenhum nulo ficou entre 0,04 e 0,06 em todos os casos: todos ficam liberais quando o alcance se aproxima do tamanho da paisagem. O variograma foi o melhor (média 0,062, máximo 0,09), seguido do paramétrico (0,066, máximo 0,10); MSR falha na região em L (até 0,23) e o toro na grade com alcance 5 ou mais (0,10 a 0,14). Resultados em `vignettes/articles/whole-tree.Rmd`. (ii) feito em 2026-11-20 por decisão do usuário, apesar do critério de 0,04 a 0,06 não ter sido atingido: `env_surrogates()` e `null = "variogram"`/`"parametric"` em `pcps_test()`, `node_test()` e `phylo_beta_test()`, com o aviso sobre alcances grandes; o estudo foi rodado de novo com o código do pacote. (iii) feito em 2026-11-20: `env_surrogates()` aceita uma paisagem (devolve paisagens com `env` ou `env_axes` trocados, com as distâncias da paisagem, inclusive toro) e matrizes `env_time` (o paramétrico move a série inteira de cada sítio). Item 1 concluído. Nulos com borda aplicados a dados reais em 2026-11-20 (`data-raw/real_data_nulls.R`, dafnídeos e anuros, com `env_surrogates(lonlat, strata)`), fechando a pendência da v0.5.2. O 2 é o mais curto e fecha uma validação antiga; o 1 é um
estudo de validade com resultado incerto e pode virar artigo próprio se nenhum nulo passar.

**Bloco B, especiação: comparação e experimentos.** Fase 4 (gen3sis, GPL, e MESS, sem
licença, só em `data-raw`, com medida de custo) e Fase 5 (persistência contra especiação nas
pontas com verdade conhecida; neutralidade regional contra taxas por ramo). Decisões
pendentes: nome final (recomendado manter `simulate_evolving_metacommunity()`; o nome
provisório `simulate_diversification()` não foi exportado) e o mapa de aninhamento como
artigo. Depois da Fase 4, o preprint pode ganhar a comparação. Da quinta limitação do
manuscrito: repetir as três replicações (Leibold et al. 2019, Urban & De Meester 2009,
Hajian-Forooshani & Chase 2026) com mais réplicas e uma varredura dos parâmetros centrais de
cada uma, em `data-raw`, para mostrar que as diferenças relatadas não dependem do ponto
escolhido. Serve de resposta pronta se a revisão pedir. Drivers em `data-raw/blockB/`
(PR #117, d9deb7f), com execução em blocos retomável; rodando no Mac a partir de 2026-11-19,
depois do lote Stan, de uma cópia congelada do pacote (custo estimado em cerca de 57 CPU-h,
3,5 a 5 h com 16 núcleos). Rodou em 7 h 20 min em 2026-11-19; resumo em
`data-raw/blockB/README.md` (`summarise_blockB.R`, tabelas em `blockB_summary.rds`). Nada
inverte uma conclusão. Mudanças a levar aos artigos de replicação e à revisão do manuscrito:
Leibold (persistência com 500 réplicas, robustez à extinção de manchas, `Vm` alto dá
neutralidade regional sem dominância local); Urban (vantagem do primeiro a chegar 76% e 74%,
não 84% e 72%; a inversão com variação permanente e a vantagem de quem evolui por mutação
dependem do contraste ambiental, o que qualifica a frase do manuscrito sobre variação
permanente); Hajian (sem mudança a 15%, não a 15% a 20%; comparação com o controle; direção
dada pela dispersão robusta à duração e à intensidade do pulso). Os três artigos de
replicação foram atualizados em 2026-11-19 (PR #122, 92a29d6), com uma seção de réplicas e
varreduras em cada um. Falta só levar ao manuscrito, na revisão do JBI: a frase sobre variação
permanente (Urban), a persistência com 500 réplicas (Leibold) e a extensão sem mudança (Hajian).

**Bloco C, inferência por ABC: robustez.** Estudos em `data-raw`, sem mudança de API, exceto
o primeiro e o oitavo: (i) estatística de sinal dos traits em `summary_stats()`, para recuperar `lambda`
(v0.3.0); (ii) PPC com modelo mal especificado de fato: filtro urbano conservado na filogenia
nos anuros, contra competição, que também testa a falha do MPD na `vignette("abc")`;
(iii) colonização limitada a partir do pool regional; (iv) validação cruzada por blocos numa
paisagem maior e sem toro; (v) árvore e paisagem variando entre simulações e detecção
imperfeita; (vi) likelihood sintética por sítio para PSIS-LOO (baixa prioridade); (vii) DOI de download do GBIF para os pools de
`anurans` (exige a conta do usuário); (viii) da quarta limitação do manuscrito: estatísticas que tornem `ne_ratio` e G
identificáveis sem ótimos ancestrais conhecidos, como a variância fenotípica dentro das
populações e a divergência entre populações relativa à distância (tipo QST), testadas com
`abc_recovery()` nos desenhos do artigo, e um guia curto de desenho amostral na `vignette("abc")`. Prioridade alta entre os itens do bloco.

**Bloco D, indivíduos e traits (ponte Leibold–Januario).** (i) vários traits por indivíduo e
traits de indivíduos que afetem a dinâmica (hoje são sorteados depois); (ii) cruzar com
dispersão e extinção de manchas (cenário 9); (iii) nulo de T_PC.PR levemente liberal na cauda
inferior (0,068 contra 0,04): sortear médias de população com o denominador observado;
(iv) `food_web_rates()` com várias espécies por nível; (v) cenário 8 de `history_rates()`
(mudanças de taxa por Poisson em eventos de reembaralhamento); (vi) `alpha` matricial no OU,
como `mvOU()`, e taxas por ramo vezes uma matriz.

**Bloco E, DAMOCLES.** Segue o plano do artigo (`dev/plano-artigo-damocles.md`): S3 com o
ClaDS2 em andamento. Variantes ainda não feitas: colonização dependente do tempo (`gamma_td`),
modelos com traits discretos (`model = 1, 2` do pacote), amostragem incompleta do pool, taxas
por época, ligação com a saída de `simulate_metacommunity()`. A suspeita de que o teste contra
partições aleatórias fosse liberal foi resolvida no S1 (calibrado com o qui-quadrado e o
bootstrap; a partição ficou só como diagnóstico).

**Bloco F, fluxo biogeográfico de Benício et al. (2021).** É a lista da antiga seção v0.8.0,
que saiu com outro conteúdo (eco-evolução). Já há partes: `observe_detection()` (detecção
imperfeita), `region_pool()`, `arrival_ages()` e `as_herodotools()`. Faltam: árvore e regiões
geradas juntas (GeoSSE), filtragem por traits discretos na comunidade, duas regiões com pools
próprios, e a camada de análise (iNEXT, Rao nos nós, GeoHiSSE, FiSSE e ES-sim, LTT por região,
DEC, MPD e MNTD em três extensões, clima e topografia, RPANDA, ClaSSE, SBEARS), tudo pelos
pacotes de referência em `Suggests`. Sobrepõe-se ao artigo do DAMOCLES e à ideia do Januário
(persistência por traits); decidir o escopo depois do S3.

**Bloco H, arquitetura genética (primeira limitação do manuscrito).** O motor assume
fenótipos Gaussianos e G constante ou resumida numa variância por população, o que parece ser
o que faltou nas replicações de modelos individuais. Plano em duas etapas: (i) G dinâmica com
equilíbrio mutação, seleção e deriva por população (redução de Bulmer transitória, não só a
permanente que já existe), que cabe no motor atual; (ii) só se (i) não bastar para a
dominância local de Leibold et al. (2019), um módulo opcional de poucos locos de efeito grande,
validado contra a replicação com indivíduos. Reprodução sexuada e ligação ficam fora. Estudo
primeiro, API depois; não bloqueia o CRAN.
Etapa (i) feita em 2026-11-22: `variance_selection = "sexual"` (desequilíbrio de ligação que a
recombinação reduz à metade por geração; `TRUE` = `"clonal"`, permanente). Nas três formas de
recrutamento do artigo do Leibold (`data-raw/bulmer_leibold2019.R`, 60 rodadas), a redução
transitória se comportou quase como sem erosão (G 45--112, d de -4 a -26; 97--100% de manchas
compartilhadas a partir de G = 200, 92--95% com a espécie pré-adaptada a partir de zero) e não
gerou dominância local. Pela regra do bloco, a etapa (ii), poucos locos de efeito grande, é a
próxima, a decidir com o usuário.

**Bloco I, metacomunidades multitróficas (pedido em 2026-11-19).** Estender o modelo de
Thompson et al. (2020) de `simulate_metacommunity()`, hoje de um nível trófico (nicho abiótico
e competição de Beverton-Holt), para vários níveis, no quadro de Guzman et al. (Ecology
Letters, doi:10.1111/ele.13162), em que os níveis diferem em dispersão e escala espacial.
Plano: (i) ler Guzman et al. e os modelos espaciais de teias já publicados antes de escrever
código, entre eles Pillai et al. (2011, PNAS, doi:10.1073/pnas.1106235108), Gravel et al.
(2011, PLoS ONE, doi:10.1371/journal.pone.0019374) e Haegeman & Loreau (2014, Ecology
Letters, doi:10.1111/ele.12214), e escolher um como alvo de replicação; (ii) matriz de
interações com sinais (consumo e competição) e resposta funcional de tipo II por par
consumidor-recurso, mantendo o nicho abiótico por espécie; (iii) dispersão, kernel e
extinção por nível trófico, para os cenários de Guzman et al.; (iv) ligação com o que já
existe: `food_web_dynamics()` (três níveis de McPeek 2017, sem espaço, eco-evolutivo) e
`food_web_rates()`, e com o motor evolutivo, para traits que mudam a interação; (v)
validação contra o caso de um nível (deve reproduzir `simulate_metacommunity()` exatamente)
e contra o modelo escolhido em (i); (vi) estatísticas por nível em `summary_stats()` e nos
artigos (alfa, beta e gama por nível, persistência). Mexe na API do motor, então pesa na
decisão de fazer antes ou depois do CRAN; recomendação: estudo e protótipo antes, API
depois do bloco A. Combina com o item 7 do bloco A (redes fluviais, em que consumidores e
recursos dispersam de modos diferentes ao longo da corrente). A Melissa Guzman é coautora de
Thompson et al. (2020) e primeira autora do quadro multitrófico; o e-mail a ela segue
pendente de decisão do usuário. Referências conferidas no OpenAlex em 2026-11-19 (só os
metadados); ler os artigos antes de citar resultados.

**Bloco J, artigo-vitrine de desenho de estudos de simulação (pedido em 2026-11-19).** Com
as adições recentes (ambiente no espaço e no tempo, nicho em vários eixos, competição por
traits, especiação, genética quantitativa, nulos com borda), falta um artigo longo no site que
mostre como desenhar um estudo de simulação a partir do tipo de pergunta, servindo de vitrine
do que o pacote faz. Plano: (i) organizar por pergunta, não por função, por exemplo: um método
detecta o processo X (erro tipo I e poder, com `run_scenarios()` e nulos)?; como a dispersão e
o nicho moldam alfa, beta e gama?; o que persiste sob ambiente que muda no tempo?; a evolução
de traits ou a especiação altera o padrão filogenético esperado?; quais parâmetros os dados
permitem estimar (ABC, `abc_recovery()`)?; (ii) em cada caso, o mesmo roteiro: pergunta,
fatores e níveis, réplicas, estatísticas resumo, nulo ou referência, custo computacional e uma
figura; (iii) apontar para os artigos que aprofundam cada caso (avaliação de métodos,
replicações, genética quantitativa, mapa de aninhamento) em vez de repeti-los; (iv) resultados
pré-computados em `data-raw/`, para o pkgdown não rodar simulações longas; (v) entrar na barra
de navegação como porta de entrada da seção "Simulating data". Depende de a API do bloco A
estar estável; recomendação: escrever depois do bloco A e antes da 1.0.0, e usar o texto
também no preprint.

**Bloco K, projeto FUNDECT/CNPq "Reconciling phenotypic evolution with metacommunity theory"
(pedido em 2026-11-21).** Objetivos em
https://provetelab.org/projects.html#reconciling-phenotypic-evolution-with-metacommunity-theory:
(1) urbanização e paisagem adaptativa de traits de várias dimensões do nicho; (2) contribuição
relativa de ecologia e evolução para a divergência de traits; (3) variação intra e
interespecífica e distribuição das espécies em gradientes de urbanização, em várias paisagens.
O que já existe: `trait_flex()` (Lepš et al. 2011: turnover, variação intraespecífica e
covariação da média ponderada da comunidade entre sítios, com regressão no ambiente),
`trait_gradient()`, `population_shifts()`, `nested_variance()`, `maladaptation()`,
`persistence()`, `dispersal_site`, `kernels` e o dataset `anurans` (cobertura urbana, duas
cidades). Itens:
(i) partição eco-evolutiva no tempo: decompor a mudança da média da comunidade entre passos de
uma simulação em sorting entre espécies, mudança dentro das espécies e covariância, e, como a
simulação guarda as médias genéticas, separar a parte genética da plástica. `trait_flex()` faz
a versão entre sítios; falta a temporal sobre `history`. Referência conferida no PDF: Ellner,
Geber & Hairston (2011, Ecology Letters 14: 603-614, doi:10.1111/j.1461-0248.2011.01616.x),
equação genótipo-fenótipo-ambiente (eqns 5 e 6): combina a abordagem de Hairston et al. (2005)
com a equação de Price em tempo contínuo e divide a mudança de uma resposta ecológica X em
evolução (Cov(m, z), mudança nas frequências de genótipos), mudança não herdável do trait
(E[Δz], plasticidade ou estrutura etária) e ambiente; estende a interações entre espécies.
Assume a "Population Mean Assumption" (X depende do trait só pela média). Na simulação, a
média genética e a fenotípica são conhecidas, então dá para calcular os três termos e
comparar com a versão retrospectiva (eqn 6), a que se aplica a dados de campo. Ainda não
verificados: Hairston et al. (2005), Price (1970, 1972) e Govaert et al. (2016). Objetivo 2.
(ii) artigo de gradiente de urbanização: sítios urbanos com menor variância ambiental e menor
conectividade (`dispersal_site`, `kernels`), com `evolution`, medindo maladaptação,
persistência e variância de traits ao longo do gradiente, para testar as hipóteses do projeto
(homogeneização e restrição da variação de traits; isolamento, menos fluxo gênico, fixação mais
rápida, maladaptação e extinção local). Ancorar só nos dados públicos de `anurans`; resultados
de manuscritos não publicados ficam fora. Objetivos 1 e 3.
(iii) o item (i) do Bloco D (traits de indivíduos que afetam a dinâmica) serve ao objetivo 3;
a análise intra e interespecífica em várias paisagens já é coberta por `trait_flex()` (com
`city` em `env`), `trait_gradient()` e `nested_variance()`.
(iv) opcional: estudo de recuperação de gradientes de seleção de Lande & Arnold (1983)
[referência a conferir] a partir de indivíduos simulados, contra a paisagem adaptativa
simulada (ótimos e `sigma`). Depende do item (i) do Bloco D, porque hoje os indivíduos são
sorteados depois da simulação e não têm aptidão própria; no lado empírico, exige uma medida de
aptidão junto da morfometria. Objetivo 1.

**Limitações do manuscrito fora do escopo.** A segunda (estrutura de estágios, restrições
metabólicas, paisagens que mudam em tempo geológico) não entra: o pacote fica no tempo
ecológico e em paisagens estáticas, e o gen3sis e o BioGEEM cobrem o tempo profundo (ver
"Visão" e o posicionamento modular).

**Bloco G, validação do `pcps.curve()`.** Estudo em `data-raw/validation_pcps_curve.R`, como
descrito na seção de métodos; só entra no pacote se tiver poder.

**v1.0.0 (CRAN).** Depois do bloco A e da Fase 4, que são os que ainda mexem na API. Os blocos
C a H são estudos e artigos que podem seguir depois do CRAN, sem mudar assinaturas. O bloco I
muda o motor; se entrar antes do CRAN, adia a 1.0.0, e depois dele pode vir como função nova
sem quebrar `simulate_metacommunity()`. O bloco J (artigo-vitrine) entra antes da 1.0.0, depois do bloco A.

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
| Traits | `trait_ou()` (um ou vários ótimos) e `paint_regimes()`, no lugar de `traits_bm()` |
| Identidade | hex sticker (opção F1) e `data-raw/hex_sticker.R` |
| Site | pkgdown em https://provetelab.org/comdynamite/, com build automático no `main` |
| Simulador estático | `simulate_static()`, na versão de desenvolvimento |
| Vignettes | `comdynamite`, `phylogenetic-structure` e `evaluating-methods` |
| Documentação | `\value` e exemplos em todas as funções exportadas; `Description` atualizada |
| Desde a 0.6.1 | `pcps_test()`, toro retangular, `summary_stats()`, `run_scenarios()`, `equilibrium_check()`, saída longa, paralelismo opcional, nulo restrito por conjunto em `trait_partition()`, piloto de identificabilidade |

## v0.1.0: primeira versão citável (concluída em 2026-11-05)

CI em macOS, Windows e Linux; 100% de cobertura de linhas; `inst/CITATION` e
`CITATION.cff`; release no GitHub com DOI do Zenodo.

- [x] Submissão ao r-universe: `comdynamite` incluído em `packages.json`, no repositório
      `diogoprov/diogoprov.r-universe.dev`.
- [x] Build confirmado em https://diogoprov.r-universe.dev/comdynamite (0.6.1, OK em
      Linux, macOS, Windows e wasm); badge e instalação pelo r-universe no README.

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
- [x] `pcps_test()`: PCPS de Duarte et al. (2016), com teste de toda a dissimilaridade e de
      cada eixo e `site_permutations`, que o pacote PCPS não tem. `torus_permutations()`
      aceita grade retangular, `symmetries` e `n_perm`. Estatísticas iguais às do PCPS 1.0.9
      (diferença 8e-15 em 20 comunidades) e p compatíveis. Em `data-raw/validation_pcps.R`,
      grade 10 x 10 sem relação real: livre rejeita 29-80% (`p_max` 3-24%), toro 4-7% com
      campos periódicos (`p_max` 1-4%) e 6-25% com bordas (`p_max` 1-10%). O PR ao PCPS
      fica para depois do aceite do manuscrito.
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
- [x] Heatmaps do artigo para visualizar os eixos de `scores`: `plot_phylo_beta()`, com `axis`.
- [x] `plot()` de `pcps_test()` com `clades` (pedido em 2026-11-10): os k clados mais
      profundos ou grupos dados por nome (famílias), em estrelas com envoltório e nome, como na
      Fig. 1 de Debastiani & Duarte (2014).
- [ ] Validação do `pcps.curve()` do PCPS (pedido em 2026-11-10; não implementado). Debastiani &
      Duarte (2014) o descrevem com a linha 1:1 da curva PSR (Diniz-Filho et al. 2012) e um nulo
      que embaralha pontas, sem validação. Duarte et al. (2018, *Ecology* 99:385–398) o
      reformularam como curva NTDR, com nulo de traits simulados por BM na árvore real e
      composição fixa, mais a curva com filogenia em estrela. Na simulação deles, a curva
      observada nunca saiu do envelope do BM, com K médio de 0,02 a 1,57, e os próprios autores
      concluem que ela não mede o sinal do trait; ela só indica quando corrigir pela filogenia a
      regressão do CWM contra o ambiente. Ressalvas nossas (não testadas): o R² sobe
      mecanicamente com o número de eixos até 1 em n − 1, e o envelope é comparado eixo a eixo,
      sem teste global (faltaria uma estatística única, como a área entre as curvas). Estudo
      proposto, com verdade conhecida em `data-raw/validation_pcps_curve.R`: variar o sinal no
      trait (λ ou o ρ de Grafen) e na distribuição (λ dos ótimos em `species_traits()`) e medir
      erro tipo I e poder da curva com o nulo de 2014 e com o de 2018, e de uma estatística de
      área; comparar com a correção filogenética do CWM de Duarte et al. (2018) e com os testes
      de CWM de Peres-Neto et al. (2017) e Zelený (2018). Só implementar no pacote se tiver
      poder; senão, o resultado vai para o artigo `whole-tree`.
- [x] Registrado em `?phylo_beta`: a matriz de MPD não é sempre euclidiana, ao contrário do
      que o artigo afirma, e o pseudo-F não exige isso. A nota antiga de que o picante não traz
      as funções que o artigo anuncia não diz quais são: fica fora do texto até ser conferida.

## v0.3.0: fluxo de estudos de simulação (concluída em 2026-11-07)

Objetivo: reduzir o código que o usuário repete. As vignettes precisaram de laços e
funções escritos à mão, o que indica o que falta no pacote.

- [x] `run_scenarios()`: grades de cenários com réplicas pareadas por árvore e
      paisagem e semente por célula, com paralelismo por `future.apply` em `Suggests`.
- [x] `equilibrium_check()`: diagnóstico de equilíbrio. Com competição filogenética e dispersão alta, a
      exclusão regional continuou até 3000 passos, e hoje nada avisa o usuário.
- [x] Saída em formato longo a partir de `history`: `as.data.frame()`, e `as_tibble()` por cima.
- [x] Benchmark de `mc_dynamics()`, em `data-raw/benchmark_dynamics.R`: 0,7 s para 50
      espécies, 100 sítios e 2200 passos, ou cerca de 19 h por 100 mil simulações em um núcleo.
      O sorteio de números aleatórios (`rmultinom`, `rpois`, `rbinom`) toma cerca de 60% do
      tempo e já é C, então código compilado não traria ganho que justifique. O ABC paraleliza
      por simulação.
- [x] `summary_stats()`: as mesmas estatísticas para dados simulados e observados
      (alfa, beta, gama, ocupação, SAD, equitabilidade, SES MPD e MNTD, decaimento com
      a distância, relação entre traits e ambiente) e, com mais de um tempo, variação
      temporal, como o CV das abundâncias. Em Guzman et al. (2022), estatísticas
      temporais aumentaram o poder explicativo em até 59%, e a dispersão dependeu delas.
- [x] Piloto de identificabilidade, em `data-raw/pilot_identifiability.R`: 30 espécies, grade
      5 x 5, árvore e paisagem fixas, sem erro de detecção; critério combinado antes de rodar,
      R2 >= 0,5 e cobertura de 85% a 95%. Com 5000 simulações, nenhum parâmetro passou: R2 de
      0,8 a 0,98 em dispersão, competição e nicho, mas cobertura de 40% a 90%, que sobe com a
      tolerância. Na sensibilidade fora do plano, com 50 mil simulações e 500 aceitas (tolerância
      0,01), os três passam com um levantamento e com série: R2 0,98 e 0,99 (dispersão), 0,88 e
      0,90 (competição), 0,91 e 0,93 (nicho), cobertura 88% a 92%. A série só estreita o
      intervalo, de 0,13 para 0,10 da amplitude na dispersão. `lambda` não é recuperado
      (R2 <= 0,07). Decisão para o ABC: seguir, com pelo menos 500 aceitas, e dar a `lambda`
      uma estatística de sinal dos traits, que `summary_stats()` ainda não tem. Falta testar com
      árvore e paisagem variando, detecção imperfeita e modelo errado.

**Concluído quando:** as vignettes puderem ser reescritas com essas funções, com bem
menos código. A `phylogenetic-structure` mostra o desenho com `run_scenarios()`; os
resultados salvos não foram refeitos. As vignettes de avaliação de nulos usam laços de outro
tipo e não mudaram.

## v0.4.0: base do hub (concluída em 2026-11-07)

Objetivo: uma interface comum, para que simulação e métodos troquem dados entre si e
com outros pacotes. As versões seguintes dependem desta.

- [x] Estrutura de dados comum: `mc_data()` e `as_mc_data()`, com comunidade, árvore, traits,
      ambiente, coordenadas, áreas e série temporal, e checagem de nomes. As análises aceitam o
      objeto no primeiro argumento, e as simulações guardam a árvore.
- [x] Árvore e distribuições como entrada: `pool` (onde cada espécie pode ocorrer) e
      `initial` (estado inicial) em `simulate_metacommunity()` e `mc_dynamics()`; a árvore
      já entrava por `species_traits()`. Cobre gen3sis, modelos de SSE e reconstruções
      biogeográficas, sem leitor próprio: o formato de saída do gen3sis não foi conferido.
      Pool que muda no tempo (expansões de área) fica para quando houver uso.
- [x] `species_traits()` aceita traits vindos de fora (`optima`, `sigma`, `r_max`), por
      exemplo do mvMORPH; os demais traits saem iguais aos sorteados sem eles.
- [x] Registrar a história dos traits ao longo dos ramos: valores nos nós em `trait_bm()`,
      `species_traits()` (`nodes = TRUE`, sem mudar os outros valores com a mesma semente) e
      `trait_ou()`; `trait_mk()` para trait discreto (Mk), no formato simmap do phytools
      (conferido contra `phytools::sim.history()`); `tip_history_metrics()` com as métricas
      de Luza et al. (2021) como o código de `calc_tip_based_trait_evo()` as calcula (nós,
      não mudanças dentro dos ramos), mais `changes` e `time_since_change` em tempo contínuo.
      O Herodotools não foi instalado (depende do daee, fora do CRAN): a equivalência vem da
      leitura do código, não de uma comparação rodada. Taxas variando entre ramos (relaxed
      BM, BAMM) como modelo gerador ficam para depois.
- [x] Integração com o Herodotools (Nakamura et al. 2024): `as_herodotools()` (entradas com os
      nomes dos argumentos dele e as áreas ancestrais verdadeiras), `region_pool()` e
      `arrival_ages()`, artigo "Known biogeographic histories for Herodotools" e
      `data-raw/validation_herodotools.R` (Herodotools 2.0.1, instalado à parte: o rexpokit do
      BioGeoBEARS não compila neste macOS, usado o binário do CRAN). Resultados: com a história
      verdadeira no lugar dos mapas, as métricas por ponta são idênticas às de
      `tip_history_metrics()` em 30 de 30 histórias. Das pontas, a correlação de Spearman com a
      verdade cai com a taxa: estase 0,97 a 0,69, transições 0,68 a 0,38; as transições nos nós
      perdem as mudanças dentro dos ramos (1,81 contra 3,01 por ponta na taxa mais alta). Idade
      de chegada com áreas verdadeiras: 0,70 a 0,89; com áreas por ML (`ape::ace`): 0,24 a 0,58,
      porque um nó profundo errado muda a idade de um clado inteiro. Herodotools fora do
      `Suggests` por não estar no CRAN.

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
- [x] `node_test()` roda nós em paralelo com `cores`. Medido em um núcleo, com 100 sítios e
      999 permutações: 0,5 min com 100 espécies, 1,7 com 300 e 3,9 com 600, linear nos nós.
      Não é lento, então `node_fit()` não foi vetorizado.
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
- [x] Estilo tidyverse na entrada e na saída, com álgebra de matrizes por dentro: tibbles
      numéricos entram, `as.data.frame()` e `as_tibble()` saem. Permutações em paralelo com
      `future.apply`, em `Suggests`, sob o `plan()` do usuário; `node_test()` segue com `cores`.
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
      dados reais. Feito em 2026-11-20 (bloco A, `data-raw/real_data_nulls.R`).
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
- [x] Nulo restrito para `c` em `trait_partition()`: `site_permutations` move `env` e `space`
      juntos, e mover MEMs de uma grade quase não muda o que explicam; no holdout o toro foi
      liberal para `c` com `env` forte. Candidato: mover só `env` nos testes de `a` e manter
      `c` com permutação livre, ou documentar `c` como descritivo com nulo restrito.
      Resultado (`data-raw/validation_fractions_null.R`): cada conjunto passa a se mover
      sozinho (`env` para `a`, `space` para `c`), o que segura o nível nos dois geradores
      (`a` .037 a .071 e `c` .030 a .080, contra .20 e .24 da permutação livre) e derruba a
      premissa de degenerescência (os MEMs movidos têm em média 30% da variância fora do span
      dos originais), mas o teste restrito de `c` só tem poder com campo no span de `space`
      (1.0; .25 a .28 com campo exponencial que `space` segura em parte, contra .98 livre),
      porque pergunta se esses mapas explicam mais do que os mesmos mapas em outro lugar da
      grade, não se o resíduo tem estrutura espacial.

## v0.6.3: base do hub (concluída em 2026-11-07)

Publica o marco v0.4.0: `mc_data()` já saíra na 0.6.2; entram traits vindos de fora em
`species_traits()`, `pool` e `initial` nas simulações, valores nos nós, `trait_mk()`,
`tip_history_metrics()` e a integração com o Herodotools (`as_herodotools()`,
`region_pool()`, `arrival_ages()`), além de `trait_k` em `summary_stats()` (marco do ABC). DOI da versão:
`10.5281/zenodo.24046171`.

## v0.6.2: fluxo de simulação, PCPS e estrutura comum (concluída em 2026-11-07)

Publica o marco v0.3.0 (`run_scenarios()`, `summary_stats()`, `equilibrium_check()`, saída
longa, benchmark e piloto de identificabilidade), `pcps_test()` com o toro retangular, as
pendências das etapas anteriores, o nulo restrito por conjunto em `trait_partition()` e o
primeiro item da v0.4.0, `mc_data()`. DOI da versão: `10.5281/zenodo.24044811`.

## v0.6.1: nulo rígido e testes das frações (concluída em 2026-11-07)

`rigid_permutations()`, enumeração fora de grade, critério para `gap`, testes das frações
`a`, `b` e `c` de `trait_partition()`, `plot_phylo_beta()`, métodos de resultados e os nomes
padronizados. Substitui a `v0.6.0` para citação: a tag `v0.6.0` traz as taxas antigas de
`patch_permutations()` e não teve release. O nulo restrito para `c` fechou depois, na
versão em desenvolvimento (item da v0.6.0).

## Inferência por ABC (concluída em 2026-11-12)

Publicada como 0.8.0 em 2026-11-12 (PR 84, tag `v0.8.0`), ao fechar o critério de conclusão com
`vignette("abc")`. Os itens em aberto abaixo seguem como trabalho futuro do bloco.

Objetivo: estimar parâmetros e comparar modelos a partir de dados observados. Depende do
piloto da v0.3.0, já feito, e da estrutura comum da v0.4.0.

- [x] Estatística de sinal filogenético dos traits em `summary_stats()`: `trait_k`, K de
      Blomberg, igual ao `picante::Kcalc()`. Piloto repetido (`pilot_identifiability_k*.rds`):
      com 5000 simulações, nada passa de novo (cobertura baixa); com 50 mil e tolerância 0,01,
      `lambda` passa por pouco, R2 0,53 e cobertura 90-91%, com intervalo de 61-63% da amplitude.
      Os outros três seguem como antes. Ressalva: o K é dos ótimos, que não se observam no
      campo; com um trait medido, a recuperação de `lambda` deve ser pior. Com 30 espécies o K
      é ruidoso; árvores maiores devem ajudar (não testado).

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

- [x] Tabela de referência: `abc_table()` (2026-11-11) sorteia parâmetros de priors uniformes
      e simula com a árvore, o ambiente e as coordenadas observados.
- [x] Detecção imperfeita na tabela de referência: `observe_detection()` (2026-11-11), com
      contagens binomiais por indivíduo e visita (Royle & Nichols 2003; Royle 2004). Em
      `data-raw/abc_detection.R` (desenho da recuperação, 20 mil simulações, três visitas,
      detecção log-uniforme entre 0,02 e 0,5, maior contagem; 9 min de simulação e 3 de
      florestas): com a detecção na tabela, R2 de 0,97 (`log_dispersal`), 0,84 (`inter_max`),
      0,91 (`sigma_mean`) e 0,45 (`lambda`), contra 0,985; 0,89; 0,92 e 0,55 com detecção
      perfeita, e `lambda` cai abaixo do critério. A própria detecção foi recuperada (R2 0,98),
      sobretudo pela abundância total; ressalva: no modelo a capacidade de suporte é fixa, então
      a abundância registrada informa a detecção diretamente; com capacidade desconhecida as
      duas se confundem (não testado). Ignorando a detecção (tabela de detecção perfeita
      aplicada aos mesmos dados), `sigma_mean` ficou enviesado em -0,36 da amplitude do prior,
      com R2 negativo e cobertura de 0,60, e `log_dispersal` caiu para R2 0,49 com viés de
      -0,11; `inter_max` e `lambda` pouco mudaram. As coberturas com a detecção na tabela
      ficaram de 0,93 a 0,985, conservadoras como antes.
- [x] Capacidade de suporte como parâmetro, para saber se a detecção continua identificável
      quando a escala da abundância é desconhecida (2026-11-11). Em
      `data-raw/abc_capacity.R` (como `abc_detection.R`, mais `log_capacity` entre -1 e 1, que
      divide os coeficientes de competição e multiplica as abundâncias de equilíbrio por 0,1 a
      10; 20 mil simulações; 8 min de simulação e 1 de florestas): o produto detecção ×
      capacidade foi recuperado com R2 0,985 e intervalo de 0,12 da amplitude, mas cada um
      separadamente só em parte: detecção com R2 0,60 (0,98 com capacidade fixa) e intervalo
      de 0,66 da amplitude, capacidade com R2 0,79, e correlação de -0,92 entre os erros das
      duas. Ambas ainda passam no critério (coberturas de 0,92 e 0,915). Os outros parâmetros
      perderam um pouco: R2 de 0,95 (`log_dispersal`), 0,76 (`inter_max`), 0,81
      (`sigma_mean`) e 0,57 (`lambda`). Conclusão: o que separa as duas, em parte, parece vir
      das estatísticas temporais e da equitabilidade (`turnover`, `cv_site`, `evenness` entre
      as mais importantes), já que a estocasticidade demográfica depende da escala da
      abundância; com um levantamento só, a separação deve piorar (não testado). Em dados
      reais, informação externa sobre a detecção (visitas repetidas com contagens por visita,
      ou um prior informativo) deve ajudar mais que o próprio ABC.
- [x] Estimação com `abcrf` (Raynal et al. 2019), em `Suggests`: `abc_fit()`, `predict()` e
      `abc_importance()` (2026-11-11). Random forests lidam bem com muitas estatísticas
      correlacionadas e ordenam a importância delas.
- [x] Escolha de modelos com `abcrf::abcrf()` (Pudlo et al. 2015) e erro de classificação
      out-of-bag: `abc_choose()` (2026-11-11). `abc` (Csilléry et al. 2012) só se a
      comparação com rejeição for útil. Em `data-raw/abc_model_choice.R`,
      competição filogenética contra nenhuma competição (10 mil simulações por modelo, mesmo
      desenho da recuperação): erro out-of-bag de 0,037 com a série e 0,039 só com o
      levantamento; em 100 pseudo-observações por modelo, acerto de 95% (competição) e 99–100%
      (sem competição), probabilidade posterior média de 0,95–0,98. O erro se concentra em
      competição tomada por neutra; nossa interpretação, não testada, é que vem dos sorteios
      com `inter_max` perto de 0, onde os modelos coincidem.
- [x] Custo, medido em `data-raw/benchmark_dynamics.R`: 0,7 s para 50 espécies, 100 sítios
      e 2200 passos, como em Guzman et al. (2022), ou cerca de 19 h por 100 mil simulações em
      um núcleo; a grade 5 x 5 do piloto levou 1,1 min por 5200 simulações em 6 núcleos.
- [x] Recuperação em dados simulados com parâmetros conhecidos: viés e cobertura
      (`abc_recovery()`, 2026-11-11). `data-raw/abc_recovery.R` repete o desenho do piloto
      grande (50 mil simulações, 200 pseudo-observações, 500 árvores; 102 min de simulação
      em 6 núcleos, 5 min de florestas; a tabela fica no rds). Com a série, R2 de 0,985
      (`log_dispersal`), 0,89 (`inter_max`), 0,92 (`sigma_mean`) e 0,55 (`lambda`), quase o
      mesmo da rejeição com tolerância 0,01 (0,98; 0,90; 0,93; 0,53). A cobertura do
      intervalo de 90% foi de 0,90–0,975: calibrada para `inter_max` e `lambda`, que passam
      no critério, e conservadora para `log_dispersal` (0,955) e `sigma_mean` (0,975), que
      falham só por cobrir demais, com intervalos mais largos que os da rejeição (0,16 contra
      0,11 e 0,31 contra 0,23 da amplitude). A rejeição com tolerância 0,01 cobriu 0,88–0,905.
      Conclusão: o ABCrf não estima melhor que a rejeição bem ajustada neste desenho; a
      vantagem é dispensar a escolha da tolerância, que no piloto mudava a cobertura de 0,64
      (0,002) para 0,90 (0,01), e nunca dar intervalos estreitos demais. `lambda` continua no
      limite (R2 0,55) e vem quase só do K de Blomberg dos ótimos (importância 0,54), com a
      mesma ressalva do piloto.

Avaliação do modelo, sempre em conjunto com a inferência:

- [x] Posterior predictive checks, com estatísticas que não entraram no ajuste: `abc_ppc()`
      (2026-11-11). No mesmo script, o modelo de competição ajustado ao levantamento
      e checado nas quatro estatísticas temporais quase não detectou o modelo errado:
      pseudo-observações sem competição foram marcadas em 5% (Bonferroni) e as com competição
      em 0%. O teste foi fraco por desenho: o modelo de competição contém o neutro
      (`inter_max` = 0), então não estava mal especificado. Falta um caso de especificação
      errada de fato, como dados com filtragem ambiental forte ajustados sem ela.
- [x] LOO sobre a tabela de referência: substituído por `abc_recovery()`, que deixa
      pseudo-observações fora e mede R2, viés e cobertura, e pelo erro out-of-bag de
      `abc_choose()`; `cv4abc()`, `cv4postpr()` e `gfit()` do `abc` não acrescentam ao que as
      florestas já dão (decidido em 2026-11-11).
- [x] K-fold por sítios ou blocos espaciais (`spatial_folds()`, `fold_summarise()`, `abc_cv()`,
      2026-11-11): ajustar com parte dos sítios e prever as estatísticas dos demais. Blocos
      deveriam evitar que a autocorrelação espacial infle o acerto. Em `data-raw/abc_cv.R`
      (desenho da recuperação, 10 mil simulações, 5 folds de 5 sítios, 100 conjuntos por
      caso; 11 min de simulação e 7 de florestas): com o modelo certo, cobertura do intervalo
      preditivo de 90% de 0,96 com blocos e com folds aleatórios, conservadora como na
      recuperação, e 1,7% de p < 0,05. Com o modelo errado (espécies seguindo um ambiente
      escondido e autocorrelacionado), a rejeição média foi quase igual (10,4% com blocos,
      10,2% aleatórios). Blocos detectaram mais em `evenness_regional` (17% contra 8%),
      `decay` (45% contra 40%) e `sad_sd` (11% contra 7%), e menos em `trait_env` (70%
      contra 87%), que detecta o erro já pelo ambiente medido. Conclusão: neste desenho
      a vantagem esperada dos blocos não apareceu; uma explicação possível, não testada, é
      que a grade de 25 sítios é um toro e os blocos de 5 sítios ficam vizinhos de todos os
      lados. Falta repetir com paisagem maior e sem toro antes de recomendar blocos.
- [x] Diagnóstico de Pareto k dos pesos de importância, para saber quando a posterior
      ponderada depende de poucas simulações: `abc_diagnose()` (2026-11-11), com
      `posterior::pareto_diags()` em vez de `loo::psis()`, o tamanho efetivo de amostra dos
      pesos e as estatísticas fora dos 99% centrais da tabela. Em `data-raw/abc_diagnose.R`
      (tabela da recuperação, 200 pseudo-observações do modelo certo e 100 do modelo sem
      competição; 5 min): no modelo certo, o tamanho efetivo foi de 141–245 e o k-hat passou
      do limiar em 0–4% dos casos, sem relação com o erro (Spearman de -0,11 a 0,24 para o
      k-hat e de 0,02 a 0,15 para o tamanho efetivo) nem com a cobertura. Com dados do outro
      modelo, o tamanho efetivo caiu para 75–153 e as estatísticas fora da tabela subiram de
      0,17 para 1,27 por observação. Conclusão: os pesos das florestas não são razões de
      importância, e o k-hat não serve como medida de confiança das estimativas; tamanho
      efetivo baixo e estatísticas fora da tabela servem como aviso de que a observação não
      se parece com as simulações. A cobertura zero de `inter_max` nos dados sem competição
      vem do verdadeiro valor estar no limite do prior (0), não de erro de ajuste.
- [ ] PSIS-LOO padrão (Vehtari et al. 2017) precisa de log-verossimilhança por observação,
      que o ABC não dá. Avaliar likelihood sintética (Wood 2010) por sítio ou por
      estatística, que daria esse termo de forma aproximada.

Limites a deixar claros: um bom ajuste não prova o mecanismo, e dados de um único
levantamento podem não informar a dispersão.

Estado em 2026-11-11: as funções estão prontas (PRs #77, #78, #80, #81 e #82): `abc_table()`,
`observe_detection()`, `abc_fit()`, `predict()`, `abc_importance()`, `abc_recovery()`,
`abc_choose()`, `abc_ppc()`, `spatial_folds()`, `fold_summarise()`, `abc_cv()` e
`abc_diagnose()`. Falta, para fechar o bloco:

- [x] Vignette com dados reais, do ajuste à avaliação: `vignette("abc")`, com os anuros de 40
      poças urbanas em duas cidades (`anurans`), pools regionais do GBIF e filtro urbano (PR 84).
- [ ] PPC com um modelo de fato mal especificado (por exemplo, filtragem ambiental forte
      nos dados e ajuste sem ela), já que o teste com o modelo neutro foi fraco por desenho.
- [ ] Repetir a validação cruzada por blocos numa paisagem maior e sem toro.
- [ ] Desempenho com árvore e paisagem variando entre simulações (hoje são fixas).
- [ ] Filtro urbano conservado na filogenia nos anuros: tolerância como trait evoluído na
      árvore (por exemplo, `lambda` próprio) em vez de limiares sorteados ao acaso. Na
      `vignette("abc")`, o modelo escolhido (competição filogenética, 0,994) não reproduz o MPD
      observado, que fica entre aleatório e agrupado (PPC p = 0,004 e 0,028; cobertura 0,25 na
      validação por blocos). Comparar os três modelos (competição, filtro conservado, os dois)
      e ver se a falha some; serve também como o PPC com modelo mal especificado acima.
- [ ] Colonização limitada a partir do pool regional (a maioria das espécies rara), a outra
      explicação da riqueza baixa que hoje só a competição consegue gerar.
- [ ] DOI de download do GBIF para os pools de `anurans` (exige a conta do usuário; hoje há
      só as consultas e a data de acesso em `data-raw/anurans_pool/`).

**Concluído quando:** uma vignette levar dados reais do ajuste à avaliação, com
cobertura medida em dados simulados.

## v0.7.0: evolução de traits e modelos nulos (concluída em 2026-11-10)

Publicada como 0.7.0 em 2026-11-10 (PR 69, tag `v0.7.0`), com todos os itens concluídos e
mergeados (PRs 30–67). Entre eles estão o PR B do DAMOCLES (PR 56), as opções novas de
`population_shifts()` (PR 53), o modelo Stan de `trait_gradient_bayes()` (PR 64), os efeitos
de fundador (PR 65) e o artigo de taxas por ramo (`vignettes/articles/rates-assembly.Rmd`,
PR 67), escrito como estudo de recuperação. O teste causal (montagem → taxas de evolução de
traits) fica para o bloco eco-evolutivo. DOI da versão: `10.5281/zenodo.24097553`.

Objetivo: ampliar os modelos de traits e os nulos com dinâmica histórica, chamando os
pacotes de referência em vez de reimplementá-los.

- [x] Mais opções de evolução de traits, multivariadas e com covariância entre traits,
      pelo mvMORPH (Clavel et al. 2015).
      Covariância entre traits feita em `trait_ou()` e `rate_matrix()`, conferida contra o
      mvMORPH (item "Funções que usam traits" abaixo).
- [x] Taxas que variam entre ramos como modelo gerador. Avaliação da literatura (2026-11-07):
      um gerador com a verdade conhecida, não um método de inferência; RRphylo, phytools,
      motmot, bifrost e corHMM são GPL e só entram em `Suggests` ou `data-raw`. O RRphylo não
      tem modelo discreto (aceita estados codificados como números, tratados como contínuos).
      BayesTraits e RevBayes adiados (binários externos; licença do BayesTraits não
      declarada); BAMM sem manutenção desde 2017.
      `trait_rates()`: modos `"shift"` (mudanças herdadas ou só no ramo), `"gbm"` (log da taxa
      por BM) e `"iid"` (lognormal por ramo), com os valores nos nós e as taxas por ramo.
      `trait_mk()`: `alpha`/`edge_rates` (Γ por ramo, como `fitgammaMk`) e `rate_classes`
      (classes ocultas, como corHMM). Validação em `data-raw/validation_rates.R`.
      `shift_rate` (mudanças por processo de Poisson ao longo dos ramos, como no prior do
      BAMM) e `relative = FALSE` (cada regime com taxa própria, como no BAMM e no bifrost).
      Validação com o BAMM 2.5.0 em `data-raw/validation_bamm.R` (binário local, GPL, só
      executado): Spearman 0,51 entre taxa média posterior e verdadeira por ramo; probabilidade
      marginal de mudança 0,18 nos ramos com mudança verdadeira contra 0,02 nos demais.
- [x] Vignette tutorial `branch-rates` (decidido em 2026-11-08 separar tutorial e artigo):
      os três geradores (`trait_rates()`, `history_rates()`, `turnover_rates()` e
      `food_web_rates()`) e a recuperação, pré-computada em `data-raw/vignette_rates.R`: o
      RRphylo segue a mudança realizada (rho 0,35–0,37) mas não a taxa verdadeira (0,06–0,08;
      teto 0,24–0,35); o BM com regimes do mvMORPH recupera a razão entre taxas (medianas 1,02,
      2,12 e 5,46 para 1, 2 e 5); a taxa do RRphylo nas pontas acompanha pouco a troca de
      parceiros (rho mediano 0,10).
- [x] Artigo (vignette) sobre taxas por ramo. Citar obrigatoriamente Hubert et al. (2015),
      modelos de especiação em metacomunidades, e Januario et al. (2025), persistência
      metapopulacional contra taxas de especiação nas pontas (sinal filogenético, mas sem
      relação com a especiação). Ponte: os dois tratam de taxas de diversificação; aqui as
      taxas são de evolução de traits, e a pergunta análoga é se a dinâmica da metacomunidade
      deixa rastro nas taxas estimadas por ramo. Para o papel da montagem de metacomunidades
      na evolução de traits, citar também Leibold et al. (2022), evolução e montagem de
      comunidades através de escalas espaciais, e McPeek (2017), *Evolutionary Community
      Ecology*. Notas de leitura e cenários de simulação em `dev/notas-artigo-taxas.md`.
      Incluir os resultados da validação e o limite de identificabilidade (uma realização por
      ramo).
      Feito em 2026-11-09: `vignettes/articles/rates-assembly.Rmd`, artigo do pkgdown, como
      estudo de recuperação (sinais plantados, não emergentes), com resumos em
      `data-raw/vignette_rates_assembly.R`. Resultados: identificabilidade das taxas por ramo;
      rastro fraco da troca de parceiros nas taxas entre espécies (rho 0,03–0,17); sinal claro
      entre populações; compensação nas estatísticas T com o nulo restrito à espécie;
      identificabilidade de c governada pela largura do kernel; poder do teste de Januario et
      al. (2025).
- [x] Ponte Leibold–Januario como eixo do artigo (detalhes em `dev/notas-artigo-taxas.md`).
      Leibold et al. (2022) apontam que os modelos macroecológicos não incorporam quanto a
      adaptação local e a montagem alteram especiação e extinção e, por elas, as filogenias.
      Januario et al. (2025) testaram uma versão disso, a persistência metapopulacional
      contra a taxa de especiação nas pontas, e acharam desacoplamento. Leibold et al. (2019)
      e Urban et al. (2020) dão mecanismos para esse desacoplamento: efeitos de prioridade
      mediados pela evolução tornam as espécies equivalentes na escala regional (neutralidade
      regional) e a evolução em geral amortece os padrões espaciais (85% dos estudos
      revisados por Urban et al.). Hipótese de trabalho (interpretação nossa): a dinâmica
      local deixa pouco rastro nas taxas entre espécies porque a divergência adaptativa fica
      entre populações, não entre espécies. O comdynamite pode testar isso com verdade
      conhecida: metacomunidade → história por linhagem (`history_rates()`) → taxas por ramo e
      traits (`trait_rates()`) → montagem → estimadores de taxa por ramo, com traits de
      indivíduos e populações (estatísticas T) para separar a escala intraespecífica da
      interespecífica. Cenários: frequência de extinção de manchas e dispersão (os dois eixos
      da Fig. 5 de Leibold et al. 2019), autocorrelação espacial do ambiente, ordem de
      chegada, as três taxas relativas da Fig. 1 de Leibold et al. (2022). Cuidado: Leibold et
      al. e Urban et al. falam de escala ecológica; Januario et al. de taxas de especiação, e
      nós de taxas de evolução de traits. O artigo deve dizer isso.
      Experimentos em `data-raw/bridge_experiments.R`, `bridge_narrow.R`,
      `bridge_bootstrap.R`, `januario_power.R` e `januario_essim.R` (notas em
      `dev/notas-januario.md`). Com clados completamente amostrados (eles usaram 189 espécies esparsas de uma
      árvore grande, não testado aqui), o PGLS de log DR sobre persistência rejeita no nível do seu
      tamanho (0,04–0,10) mesmo com elo causal forte, enquanto o ES-sim (Harvey & Rabosky 2018)
      sobre τ detecta 0,59, 0,89 e 0,98 para b = 0,25, 0,5 e 1, com tamanho 0,015.
- [x] `food_web_dynamics()`: modelo eco-evolutivo de três níveis tróficos do capítulo 3 de
      McPeek (2017), eq. 3.7–3.17 e 3.20–3.21, reimplementado das equações publicadas (o
      código MATLAB do autor não tem licença declarada; em novembro de 2026 estava fora do ar em
      enallagma.com e não arquivado, então a validação é contra o livro). Abundâncias e
      médias dos traits integradas juntas por RK4 de passo fixo em R base; gradientes de
      seleção por diferenças centrais, o que permite as três formas de coeficiente de ataque
      (gaussiana, logística, Michaelis-Menten) e a omnivoria sem derivadas à mão. Níveis
      ausentes e introduções em tempos dados (`introductions`), para testar a previsão de
      que a evolução dispara após uma invasão e desacelera até o equilíbrio e que o predador
      enfraquece a seleção sobre o recurso; a saída guarda os gradientes e |dz/dt| por nível,
      para virar taxa por linhagem. Validação em `data-raw/validation_food_web.R` (fig. 3.2,
      3.3 e 3.11). Duas leituras da eq. 3.20 documentadas na ajuda: o termo de perda do
      recurso para o predador sem `w`, e `m` como função dos traits do consumidor e do
      predador (box 3.1).
      Conferência equação a equação (subagente, nov. 2026): gradientes numéricos iguais aos
      analíticos da eq. 3.20 (diferença relativa ≤ 2e-7); Fig. 3.2 e painel A da Fig. 3.11
      reproduzidos quase exatamente. Discrepâncias atribuídas ao livro: legenda da Fig. 3.11
      com z̃c_R = 20 nos painéis A–D (o recurso sozinho fica em 12) e painel E (logístico)
      irreprodutível com os parâmetros da legenda (o consumidor não persiste no traço
      plotado). Corrigido δ = 0,02 nos painéis I–L na validação. Contra o código MATLAB dos
      modelos recurso–consumidor do suplemento de McPeek (2017, *Am Nat*): funções de aptidão
      e gradientes idênticos e trajetórias da Fig. 3.2 a 1e-4 do `ode45` do autor; sem
      predador nesse código. Acrescentado `trait_floor` (piso dos traços em zero, como no
      código do autor; padrão com Michaelis-Menten).
      `food_web_rates()` traduz as taxas por nível em `edge_rates` para `trait_rates()`
      (cenário 4 das notas), na janela final dos ramos terminais, como `turnover_rates()`.
      Falta: várias espécies por nível (cap. 5).
- [x] Função que traduz a história ao longo da árvore em taxas por ramo: `history_rates()`,
      com taxa por estado (`"state"`, como O'Meara et al. 2006), explosão após cada mudança
      que decai exponencialmente (`"burst"`, cenários 1 e 2 das notas) e dependência do número
      de linhagens no mesmo estado (`"sympatry"`, o DDexp com geografia de Drury et al.
      2016); `trait_rates(edge_rates = )` recebe as taxas. `"sympatry"` validado contra o
      RPANDA em `data-raw/validation_sympatry.R`: matrizes de covariância iguais às do DDexp
      com geografia (diferença máxima 3e-15) e r recuperado por `fit_t_comp()` (média -0,30
      para -0,3 verdadeiro, 0,005 para 0). `CreateGeoObject()` falha em mapas sem mudança de
      região (bug do RPANDA, objeto `old.Dist`); essas histórias são sorteadas de novo.
      Ligação com `simulate_metacommunity()` (cenário 5): `partner_turnover()` conta entradas
      e saídas de parceiros por sítio e passo, e `turnover_rates()` as converte em taxas nos
      ramos terminais, na janela `window` do fim de cada ramo (a simulação roda em tempo
      ecológico só entre as pontas). Falta o cenário 8 (mudanças de taxa por Poisson
      concentradas em eventos de reembaralhamento e compartilhadas entre linhagens simpátricas).
- [x] Indivíduos a partir das simulações (variação intraespecífica e adaptação local por sítio)
      e estatísticas T de Violle et al. (2012) via `cati::Tstats()` em `Suggests` (cenário 10
      das notas). `sample_individuals()` sorteia n indivíduos por população (ou Binomial(N,
      `fraction`)) com trait = ótimo da espécie + fração `local_adaptation` da distância ao
      ambiente do sítio + desvio da população (`population_sd`) + ruído individual
      (`intraspecific_sd`). `tstats()` calcula as três razões em R base (sem nulos, para
      estudos com muitas réplicas) ou via `cati::Tstats()` (nulos); os valores observados
      batem. O cati casa nomes de sítio por regex (`site1` pega `site10`): os sítios vão com
      rótulos de largura fixa e voltam com o nome original. `data-raw/demo_tstats.R`: a
      adaptação local sobe T_IP.IC (0,02 → 1,0) e baixa T_IC.IR e T_PC.PR (0,54 e 0,89 → <0,01).
      Isso contraria a previsão do cenário 10 (T_IP.IC cair): adaptação ao ambiente do sítio
      torna as espécies co-ocorrentes parecidas; T_IP.IC só cai com divergência entre
      co-ocorrentes (deslocamento de caráter).
      Segunda rodada: deslocamento de caráter em `sample_individuals(displacement = ,
      displacement_width = )`, equilíbrio entre seleção em direção à média local e repulsão
      por um núcleo gaussiano de competição ponderado pela abundância dos competidores; nulos
      em R base em `tstats()` (os três de Taudiere & Violle 2016, com SES e p-valor), de modo
      que o cati só serve para conferir. O nulo de T_IP.IC do cati junta permutações de
      sítios diferentes (a permutação i do sítio s usa indivíduos de outro sítio), então os
      SES de T_IP.IC do cati estão errados. Demo: o deslocamento baixa T_IP.IC e sobe T_PC.PR
      e contrabalança a adaptação local (h = 0,9: T_IP.IC 0,61 → 0,011 e T_PC.PR 0,010 → 0,59
      com `displacement = 10`). O SES de T_IP.IC já é muito negativo sem deslocamento: o nulo
      testa se as populações do sítio diferem, não se divergiram. Falta: vários traits por
      indivíduo (agora que `assemble_community()` tem vários traits), traits de indivíduos que
      afetam a dinâmica (hoje são sorteados depois) e cruzar com dispersão e extinção de
      manchas (cenário 9) para a ponte Leibold–Januario.
      Terceira rodada (identificabilidade): `tstats(null_ip = "species")` troca as médias das
      populações entre populações da mesma espécie (SES ≈ 0 sem processo, > 0 com
      convergência, < 0 com divergência; conservador com n desigual, nunca anticonservador;
      trocar populações inteiras testaria a filtragem de abundâncias e deu 22% de rejeição sob
      H0). `population_shifts()` ajusta h e c juntos (regressão dentro de espécie no ambiente
      e no empurrão previsto pelos competidores, metades sorteadas contra o erro amostral, w
      por perfil refinado com `optimize()`). Piloto ABC
      (`data-raw/pilot_tstats_identifiability.R`, metacomunidades de teste fora da referência,
      500 pseudo-observações): só as razões dão R² 0,87–0,89 para h e 0,56–0,58 para c;
      com `population_shifts()`, 0,99 e 0,95–0,96, recuperáveis na tolerância 0,05. Auditoria
      independente (Fable) antes da rodada final mudou o desenho do piloto e o refinamento de
      w. Limites: modelo do piloto é o mesmo que gera os dados, sem `population_sd`, com n = 10
      e `intraspecific_sd` = 2 fixos; `c` atenua com c alto e com `fraction`.
      Quarta rodada (substitutos do cati): auditoria completa do cati 0.99.6 achou bugs em
      `ComIndexMulti`, `ComIndex` (regional.pop, ordem de `sp`), `sum_Tstats`, rótulos de
      `Tstats`, p-valores e `ses()` (cati#13 a #21). Reescritos a partir das definições:
      `tstats_summary()`, `trait_flex()` (Lepš et al. 2011), `nested_variance()` (checa a
      ordem do aninhamento) e `community_index()` (quatro nulos, multivariado). Conferidos
      contra o cati onde ele está certo e na reanálise do capítulo 1 do Leandro.
      Auditoria Fable das funções novas: corrigido o denominador do nulo de T_IC.IR em
      `tstats()` (recalculado a partir dos valores sorteados, anticonservador; o cati faz o
      mesmo); `trait_flex()` com NA em `env`, variável chamada `y` e fórmula sem intercepto.
      Em aberto: o nulo de T_PC.PR (método publicado, denominador recalculado) tem tamanho
      um pouco acima do nominal na cauda inferior (0,068 contra 0,04); pensar num nulo que
      sorteie médias de população diretamente, com o denominador observado.
      Quinta rodada (PRs 47 e 49, ferramentas para dados de indivíduos): métodos `plot()` em
      ggplot2 para `tstats()`, `community_index()`, `trait_flex()` e `nested_variance()`
      (intervalos nulos aninhados ou SES; `null_intervals()` devolve a tabela);
      `trait_densities()` e `plot_population_traits()` (densidades por população e comunidade,
      ideia do TPD de Carmona et al. 2016, 2019; aviso quando a largura de banda única ajusta
      mal, com `highlight`); métricas clássicas de vizinhança (`mnnd`, `sdnnd`, `mst_length`,
      `trait_range` etc.); esforço amostral com `sample_individuals(cover, min_n)`,
      `subsample_individuals()` e `trait_rarefaction()`. Testadas com os dados reais do
      Leandro (983 anuros); `population_shifts()` se aplica a esses dados (espécies com ≥ 2
      populações de ≥ 2 indivíduos).
- [x] Funções que usam traits, como `assemble_community()`, com mais de um trait e com
      restrições na estrutura entre traits, como em `mvSIM()` e nos ajustes do mvMORPH:
      sem covariância (`"diagonal"`), taxas iguais (`"equal"`), variâncias ou correlações
      iguais (`"variance"`, `"correlation"`), matrizes proporcionais (`"proportional"`) ou
      compartilhadas (`"shared"`), e matriz definida pelo usuário. Definir como cada processo
      usa vários traits: distância multivariada na similaridade limitante, nicho
      multidimensional no filtro.
      Feito (estrutura entre traits): `trait_ou(sigma2 = )` aceita matriz de taxas ou lista
      nomeada com uma matriz (ou taxa) por regime, e `regimes` aceita um mapa simmap (de
      `trait_mk()`), com regimes, ótimos e taxas mudando dentro dos ramos; com `sigma2`
      escalar os valores são os mesmos de antes com a mesma semente. `rate_matrix()` monta as
      estruturas do `mvBM()`: `"default"`, `"diagonal"`, `"equaldiagonal"`, `"equal"` e, entre
      regimes, `"proportional"`, `"shared"` (mesmos autovetores), `"variance"`, `"correlation"`.
      Todas fazem sentido aqui porque o pacote já tem regimes (`paint_regimes()`, simmap). A
      matriz de índices inteiros do mvMORPH (restrições definidas pelo usuário) só serve para
      ajuste: na simulação o usuário passa a matriz. Validação em
      `data-raw/validation_mvtraits.R` (mvMORPH 1.2.3): covariâncias das pontas dentro de 2,9
      erros de Monte Carlo de `sum_k kronecker(R_k, C_k)` (mvSIM: 2,8), OU dentro de 3,1, e
      `mvBM()` recuperou as matrizes (mediana de 30 ajustes a menos de 9% de cada elemento).
      Ficou de fora: `alpha` matricial e `theta` por trait no OU (como `mvOU()`), e taxas por
      ramo de `trait_rates()` vezes uma matriz (estrutura proporcional por ramo), para não
      mexer em `R/traits_rates.R` enquanto `history_rates()` é ligada à metacomunidade.
      Feito (uso de vários traits, pedido em 2026-11-08): `assemble_community(distance = )`
      com `"euclidean"` (padrão, como antes), `"standardized"` e `"mahalanobis"` (covariância
      do pool ou `trait_cov`), na similaridade limitante e no filtro; a ajuda justifica a
      escolha (euclidiana conta duas vezes eixos correlacionados; Mahalanobis trata cada
      direção independente como um eixo; padronizada quando cada trait é um eixo distinto).
      Nicho com vários eixos ambientais em `growth_rates()` e `simulate_metacommunity(env = )`,
      com uma lista de traits por eixo e o produto das respostas, como em `simulate_static()`.
      Falta: nicho com correlação entre eixos (matriz de largura por espécie) na dinâmica;
      competição a partir da distância entre traits (hoje `competition_matrix()` é aleatória
      ou filogenética); paisagens que gerem vários eixos ambientais de uma vez.
      Auditoria (2026-11-08): corrigidos `theta` nomeado reciclado, simulação com vários eixos
      sem árvore (e `as_mc_data()`/`summary_stats()`/`sample_individuals()` a jusante), traits
      com NA na montagem, simmap em outra ordem de ramos; `theta` por trait; classes
      `rate_matrix` e `mc_assembly`; `plot_traits()`, `plot.rate_matrix()`, `plot_assembly()`;
      GAUSE 23x mais rápido com a mesma ordem; artigo "Correlated traits and assembly in trait
      space". Falta ainda: `alpha` matricial no OU (como `mvOU()`).
- [x] Modelo nulo dinâmico do DAMOCLES (Pigot & Etienne 2015): colonização e extinção
      local ao longo da árvore. Feito em `R/damocles.R`: `simulate_damocles()`,
      `damocles_loglik()`, `damocles_fit()` (ML em base R) e `damocles_test()` (MPD/MNTD contra
      DAMOCLES e contra sorteio aleatório); validado contra o pacote DAMOCLES em
      `data-raw/validation_damocles.R`. Falta: as variantes do pacote (colonização dependente
      do tempo, gamma_td; modelos com traits binários ou ternários, `model = 1, 2`;
      colonização dependente de traits ou de parentesco das simulações do artigo); amostragem
      incompleta do pool; ligar à saída de `simulate_metacommunity()` (vários sítios como
      comunidades locais da mesma árvore). Ajuste com taxas compartilhadas entre comunidades
      (`shared = TRUE`) feito em 2026-11-08; falta perfil de verossimilhança ou intervalo para
      as taxas, e taxas que variam entre comunidades com uma covariável. Regra de especiação
      relaxada (`sympatry`, probabilidade de as duas filhas ficarem presentes) feita em
      2026-11-08. Taxas por regime de ramo e com covariável por ramo (`regimes`,
      `edge_covariate`, `vary`), `clade_regimes()`, `damocles_lrt()` e `confint()` por perfil
      feitos em 2026-11-08 (PR A das ressalvas). PR B feito em 2026-11-08:
      `damocles_samples()` (lista de configurações com resumo), `damocles_partition_test()` e
      `random_regimes()` (partições aleatórias de mesma complexidade), `clade_regimes()` com
      nós, `dr_rates()` (DR de Jetz et al. 2012), conversão de BAMM e ClaDS em
      `data-raw/damocles_bamm.R` e `data-raw/damocles_clads.R`; `max_overlap` em
      `random_regimes()` (sem ele o teste contra partições aleatórias não tinha poder). Falta:
      checar se o teste contra partições aleatórias é levemente liberal (4/40 sem efeito) com
      mais réplicas.
      Próximos, em ordem de
      facilidade: taxas por época (quadratura do fator integrante, sem EDO); trait discreto coevoluindo com a presença (estado
      presença × trait); traits contínuos e colonização dependente de parentesco só por
      simulação.

- [x] Análise de gradiente de traits de Ackerly & Cornwell (2007) (pedido em 2026-11-09):
      `trait_gradient()`, `trait_gradient_weights()`, `trait_gradient_test()`,
      `trait_gradient_phylo()` e `plot()`. Reimplementada das equações (o script do material
      suplementar não tem licença; só foi lido e usado para conferir: `leave_one_out = FALSE`
      bate com ele a 1e-14). Três correções: médias de sítio sem a espécie focal; pesos de
      sítio e de espécie separados; nulo e erro (troca de médias entre espécies como nulo
      padrão, curveball só para a correlação β–α, GLS por REML com erro de bootstrap em BM/OU).
      Validação: o curveball rejeita 26–66% nos testes de variância sem triagem de traits, por
      isso deixou de ser o padrão; o nulo por espécies fica em 4,7–7,3%, mas sobe a 10–17%
      quando trait e nicho têm sinal filogenético (associação herdada). Auditoria Fable feita
      e corrigida.
- [x] TGA, extensões (pedido em 2026-11-09): nulo `"phylo"` com as permutações filogenéticas
      de Lapointe & Garland (2001), conferidas com o PDF (k ≥ 1, valor de i vai para j), e nulo
      `"bm"`; validação: com trait e nicho filogenéticos o erro cai de 14–23% para 2–11%, e o
      `"bm"` falha para traits sem estrutura filogenética. `trait_gradient_brms()`: modelo
      bivariado em brms com β e média latentes, parte filogenética e não filogenética, λ com
      intervalo, priors explícitos, diagnósticos e `pp_check()`; OU com α dado (sugestão do
      usuário: A da covariância OU com raiz fixa). Com 18 espécies λ é bimodal e as cadeias não
      misturam. Validação com ~43 espécies (20 réplicas): λ = 1 recuperado em 8/10 (5 passam
      em todos os diagnósticos); λ = 0 recuperado em 8/10, mas nenhuma passa (ESS 180–400);
      4 réplicas com cadeias divididas. Auditoria Fable corrigida (nomes com espaço, Kish,
      diagnósticos de λ, loo por resposta, β NA para espécies sempre sozinhas).
      Validação 2 (4000 iterações, pesos de ocorrência, ~43 espécies): R-hat ≤ 1,03, mas
      divergências na maioria; λ = 0 de novo perto de 0 (10/10 com divergências); λ de Pagel
      0,5 NÃO é recuperado: medianas 0,007–0,81 (mediana 0,04), intervalo contém 0,5 em 12/20 e
      exclui 0 em 3/20; só 4/20 passam. Documentado em Limits.
- [x] TGA bayesiana reparametrizada em Stan próprio (pedido em 2026-11-09):
      `trait_gradient_bayes()` substitui `trait_gradient_brms()` (nunca na main); τ total e
      λ ~ beta(1, 1) por resposta, efeitos não centrados, cmdstanr (BSD) em Suggests, modelo
      em `inst/stan/` compilado no cache do usuário. Validação (60 réplicas, ~43 espécies):
      9–30 s por ajuste, 1 divergência em 60; λ = 0 → mediana 0,07; λ = 0,5 → 0,24, intervalo
      contém 0,5 em 17/20 (brms: 12/20, mediana 0,04); λ = 1 → 0,94. Falhas de diagnóstico
      só por ESS < 400.

- [x] Rodada de UX antes de fechar a 0.7.0 (pedido em 2026-11-08; concluída em 2026-11-09):
      equações da ajuda em LaTeX pelo KaTeX (PRs 46 e 48); vignette `branch-rates` (PR 42).
      PR 52: `plot_rates()` (taxas por ramo pintadas na árvore) e `plot()` do PCPS; mensagens de
      erro que nomeiam o argumento (`check_tree()`, `check_count()`); `@family` em 79 funções.
      PR 54: renomeações com aviso via `.Deprecated()` até a 1.0.0 (`traits_ou()` →
      `trait_ou()`, `fit_damocles()` → `damocles_fit()`, `nperm` → `n_perm`, `root_value` →
      `root`, `trait_mk(root)` → `root_state`, `$edge_rate` → `$edge_rates`, `colour` →
      `colours`; `tstats_summary(level)` passa a ser a cobertura). PR 55: todos os plots em
      ggplot2 (`Suggests`), retornando objetos ggplot com as tabelas em atributos.

**Concluído quando:** cada opção tiver teste e exemplo em vignette.

## Especiação acoplada à metacomunidade (número a definir, depois da v0.7.0)

Objetivo: fechar o limite da ponte Leibold–Januario. Hoje a árvore é dada e a especiação não
depende da metacomunidade. Plano completo, com o modelo, os modos de especiação, a validação e
as fases, em `dev/plano-especiacao.md`.

- [x] Fase 1: traits por população (espécie × sítio) com seleção local e fluxo gênico em
      `mc_dynamics()`, sem especiação (2026-11-12): `evolution_control()`, com deriva,
      fundadores, plasticidade e passos por geração; testes contra Lande (1976) e o equilíbrio
      migração–seleção. Neutralidade regional reproduzida qualitativamente em
      `data-raw/evolution_leibold.R`: só com extirpação, dispersão baixa e adaptação rápida
      (sorting de 0,99 para 0,04 com G = 50, contra 0,90 no controle sem resposta).
- [x] Fase 2a (feita em 2026-11-12): `simulate_evolving_metacommunity()`, `speciation_control()`
      (`"geographic"` com regras de clusters e gradual; `"ecological"`), árvore completa e
      reconstruída, extinção global, `as_metacommunity()`; validação 3 sem apoio (ρ ≈ 0).
      Revisada antes da 2b: brotamento, ligação por taxa de migração com persistência
      (`migration_min`, `reset_generations`), modo ecológico com persistência e fluxo reduzido
      (`min_generations`, `migration_max`) e competição pela sobreposição de nicho.
- [x] Fase 2b (feita em 2026-11-12): `"point"`, `"fission"` (partição espacial ou aleatória),
      `jump_var`, conclusão prolongada (fixa ou exponencial) com agrupamento das linhagens em
      espécies boas; validações 1 (Ewens) e 2 (γ e PBD).
- [x] Fase 2: divergência entre populações e especiação (iniciação por `"point"`,
      `"fission"`, `"geographic"` ou `"ecological"`, seguindo os quadros de Hubert et al.
      2015, e uma camada de conclusão prolongada como no PBD), com registro filogenético
      (árvore completa e reconstruída, tempos de iniciação e de conclusão) e extinção global;
      validação no limite neutro, contra o PBD e contra Desjardins-Proulx & Gravel (2012a).
- [x] Fase 3a (feita em 2026-11-12): séries por linhagem, `persistence()` (ε e τ de Januario et al.
      2025), `sample_individuals(source = "evolution")` para as estatísticas T e
      `landscape_matrix()` (kernel por matriz, para os dados do Marcos).
- [x] Fase 3b (feita em 2026-11-12): registro das mudanças de presença (`ranges`) e
      `range_history()`, simmap com estado = combinação de áreas ocupadas, na árvore completa
      ou reconstruída; `as_herodotools()`, `region_pool()` e `arrival_ages()` com áreas
      combinadas; `history_rates(overlap = TRUE)` na simpatria.
- [x] Bloco próprio, depois da Fase 3: G multivariada (vários traits com covariância
      genética, B ∝ W), pedida para simular os dados do Marcos (ECY25-0086; decidido em
      2026-11-12).
  - [x] Parte 1 (sem especiação): um trait por eixo de nicho, G e E como matrizes em
        `evolution_control()`, seleção G(Ω+P)⁻¹d, deriva MVN, resposta correlacionada nos
        eixos com `sigma = Inf`; `genetic_matrices()` (proporcional, independente, PCs
        compartilhados); `sample_individuals()` com um traço por coluna. Testes: G = 0,
        deriva Var(Δz) = G Σ 1/Nₑ, resposta em um passo, razão G21/G11.
  - [x] `maladaptation()` (desvio e lag load por população; Mahalanobis com vários traits)
        e artigo "Quantitative genetics in the metacommunity" com GIFs (carga de migração,
        resgate evolutivo, orientação da G, prioridade evolutiva; Leibold et al. 2019, Urban &
        De Meester 2009, Cabral e colaboradores).
  - [x] Parte 2: especiação com vários traits em `simulate_evolving_metacommunity()`
        (especiação ecológica pela distância de Mahalanobis sob P, salto MVN com `jump_var`
        matricial, sobreposição de nicho produto dos eixos); `env_time` com vários eixos;
        `drift_test()` (Ackermann & Cheverud 2004) e validação B ∝ W em
        `data-raw/validation_drift.R` (deriva contra seleção divergente no eixo menor da G).
- [x] Mapa de aninhamento e replicações (decidido em 2026-11-12; antes da Fase 4): tabela de
      quais configurações recuperam cada modelo publicado e quais mecanismos ficam de fora
      (genética individual de Urban & De Meester 2009, distúrbio de Leibold et al. 2019,
      restrições metabólicas, estágios e dinâmica geoclimática de Cabral e colaboradores/BioGEEM,
      tempo profundo do gen3sis). Replicar, nesta ordem: (i) a convergência de nicho regional de
      Leibold et al. 2019 (o artigo de genética quantitativa só mostra a monopolização); (ii) o
      efeito do atraso de colonização de Urban & De Meester 2009. Prioridade sobre o gen3sis,
      que é bom para tempo profundo e menos para tempo ecológico. O posicionamento (populações e
      genética quantitativa em tempo ecológico, módulos testáveis com a verdade conhecida;
      Cabral, Valente & Hartig 2017) vai para o README e para o preprint.
      Artigo "Nesting map" (2026-11-13): seis modelos, status, matriz mecanismo × modelo
      (`data-raw/figure_nesting_map.R`); posicionamento no README (mesmo PR).
      Feito (2026-11-13): rascunho do mapa em `dev/nesting-map.md` (seis modelos, fontes por
      item; BioGEEM só pelo preprint, Tabela S1 de Urban & De Meester 2009 não lida);
      `disturbance` no motor; replicação (i) em `data-raw/replication_leibold2019.R` e no
      artigo "Replicating Leibold et al. (2019)": neutralidade regional reproduzida com G ≥ 100;
      dominância local mais fraca e deriva regional mais rápida que no modelo deles (hipótese:
      imigrantes misturados à média). Testado (`immigrant_selection`, `data-raw/mechanism_leibold2019.R`):
      a seleção dos imigrantes não explica nenhuma das diferenças; deriva e amostragem de
      fundador das médias explicam a maior parte da deriva regional mais rápida; a dominância
      local fraca parece vir do G que não depende do tamanho da população (a testar com G
      dependente de N). Testado (`mutation_var`, variância por população com mutação e deriva):
      nenhuma espécie se perde, a neutralidade regional surge, mas a dominância local some
      (90-99% de manchas compartilhadas com Nm = 10); a dominância local com G constante vem da
      deriva das médias. A do modelo deles fica sem explicação (candidatos não testados:
      recrutamento em loteria, gerações sobrepostas, variância a partir de mutações novas).
- [x] Três mecanismos (pedido em 2026-11-13): (i) recrutamento em loteria nos microssítios que
      vagam (`recruitment = "lottery"`, `microsites`), (ii) gerações sobrepostas (`mortality`
      < 1) e (iii) variância construída por mutação e erodida pela seleção (`variance_selection`
      em `evolution_control()`, redução de Bulmer permanente). Testados em fatoriais
      (`data-raw/mechanisms_urban2009.R`, `data-raw/mechanisms_leibold2019_recruitment.R`,
      `data-raw/figures_mechanisms.R`) e seções nos dois artigos. Urban: a erosão pela seleção
      restaura o efeito de prioridade da variação permanente; gerações sobrepostas restauram a
      vantagem do melhor dispersor com mutação; loteria com gerações discretas não muda nada.
      Leibold: nenhum dos três produz dominância local; o resultado depende do regime de
      variância (G constante: perda de espécie; mutação a partir de G: manchas compartilhadas;
      a partir de zero ou erodida: species sorting). Em aberto: V_m maior que 0,5.
- [x] Replicação (iii), Hajian-Forooshani & Chase 2026 (PR #106, 2026-11-17): `env_pulse()` e o
      artigo `replication-disturbance`, estatisticamente indistinguível das simulações
      arquivadas pelos autores; o mapa de aninhamento passou a sete modelos.
- [x] Replicação (ii), Urban & De Meester 2009 (`pool_from`, `data-raw/replication_urban2009.R`,
      artigo `replication-urban`): a espécie que evolui por mutação vence o sítio novo em
      qualquer ordem de chegada e a evolução nas duas torna o sítio mais monopolizado, mas
      sem favorecer quem chega primeiro; variação permanente só reproduz com G constante (carga
      de variância no `mutation_var`); na grade de dispersão a evolução faz o pior dispersor
      vencer (região de maladaptação por fluxo gênico), sem a região de monopolização (sem
      evolução, o melhor dispersor vence). Correção no caminho: sítio fechado por `pool_from`
      não é mais sumidouro de emigrantes (eles ficam em casa); a primeira versão atribuía o
      resultado da dispersão ao custo da emigração, por engano. A perda
      da fonte por fluxo gênico ocorre com dispersão baixa, não alta. Os três mecanismos da
      lista para depois são também os candidatos aqui. Em aberto: o mapa como artigo.
- [ ] Fase 4: comparação com o gen3sis (GPL) e o MESS (sem licença declarada), só em
      `data-raw`, e medida de custo.
- [ ] Fase 5: experimentos (persistência contra especiação nas pontas com verdade conhecida;
      neutralidade regional contra taxas por ramo) e artigo.
- [x] Interface: função nova (decidido em 2026-11-08), com `as_metacommunity()` para os
      métodos atuais; nome provisório `simulate_diversification()`.
- [ ] Decidir o nome final. Modos e número de versão já saíram (0.8.0, quatro modos).
      Recomendado em 2026-11-17: manter `simulate_evolving_metacommunity()`.
- [x] Avaliar as oito sugestões do ECY25-0086 para a Fase 1 (deriva explícita, G
      multivariada, plasticidade separada da genética, fundador mecanístico, ótimo variável no
      tempo, passos por geração, amostragem alocrônica e testes analíticos). Conferido em
      2026-11-17: sete feitas; a 5 (ótimo variável no tempo) está feita sem a saída de
      "tempo desde a mudança", que foi para o bloco A.

**Concluído quando:** as fases 1 a 4 tiverem testes, validação em `data-raw` e uma vignette.

## v0.8.0: traits, filogenia e biogeografia

Objetivo: transformar o desenho de Benício et al. (2021) em um fluxo do pacote e testar,
com a verdade conhecida, se ele separa tempo para especiação, taxas de diversificação
e conservação de nicho. Scripts e dados do artigo estão no Dryad
(https://doi.org/10.5061/dryad.sf7m0cg65).

Simular o cenário:

- [ ] Árvore e regiões geradas juntas, por exemplo com um simulador de GeoSSE, em cenários
      que variem especiação, extinção, dispersão e tempo de colonização por região.
- [ ] Filtragem por traits discretos, como os modos reprodutivos. `trait_mk()` (v0.4.0)
      já simula o trait com a história; falta o modelo de comunidade usá-lo.
- [ ] Duas regiões com pools próprios e dispersão limitada entre elas, usando a
      entrada da v0.4.0.
- [ ] Amostragem imperfeita: sítios, dias de coleta e espécies não detectadas.

Analisar, chamando os pacotes de referência sobre a estrutura comum:

- [ ] Riqueza por rarefação e extrapolação (`iNEXT`).
- [ ] Sinal filogenético em traits discretos: Rao decomposto nos nós (`adiv`, `ade4`).
- [ ] Diversificação por região: GeoHiSSE (Caetano et al. 2018), no lugar do GeoSSE
      usado no artigo, e testes baseados nas pontas, como FiSSE (Rabosky & Goldberg
      2017) e ES-sim (Harvey & Rabosky 2018).
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

## v0.9.0: mais processos (concluída em 2026-11-17)

Objetivo: ampliar os cenários do simulador. Escolher pelo uso real.

- [x] Vários eixos ambientais: `env` com uma coluna por eixo e `env_time` em lista ou array.
- [x] Paisagens irregulares ou lidas de coordenadas: `landscape_points()`.
- [x] Distúrbio: `disturbance` (perda de populações) e `env_pulse()` (pulsos no ambiente, com
      duração e, desde a PR #107, frequência via `every`).
- [x] Recrutamento em loteria, gerações sobrepostas e abertura de sítios no tempo
      (`recruitment = "lottery"`, `microsites`, `mortality`, `pool_from`).
- [x] Kernel Gaussiano (`disp_shape`, PR #107) e dispersão variável no tempo (`dispersal_time`,
      `kernels`, `kernel_step`, PR #108), com o artigo "Dispersal and disturbance through time"
      (PR #110, exemplo de pulso de inundação).

Para depois da 0.9.0 (proposto em 2026-11-17):

- [ ] Null models restritos com bordas para preditores contínuos e para o simulador. Em
      `pcps_test()`, o toro numa grade com bordas rejeitou 6% a 25% sem relação real.
      `patch_permutations()` cobre só eventos; candidato: Constant Moran permutation
      (Arthur 2024), a testar.
- [ ] Gerador de ambiente temporal fiel ao de Thompson et al. (2020), como alternativa
      a `env_over_time()`.

**Concluído quando:** os itens escolhidos tiverem testes e exemplo em vignette. Fechada
em 2026-11-17 com os itens marcados, que têm testes e artigos.

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

- **Ordenação PCPS na vignette `whole-tree`.** Pedido em 2026-11-08: além de PCPS 1 contra o
  ambiente, mostrar a ordenação das comunidades no plano PCPS 1 × PCPS 2 com o ambiente como
  tamanho (ou cor) dos círculos, em R base ou em ggplot2 (só na vignette, via `Suggests`).
  Ver se `pcps_test()` guarda os escores do segundo eixo e as posições das espécies (biplot).
- **Imigração na grade.** Decidido em 2026-11-07: argumento `immigration` em
  `simulate_metacommunity()` e `mc_dynamics()`, com `"poisson"` como padrão (código original,
  sementes antigas iguais) e `"multinomial"` conservando os dispersores, como o mcomsimr. No
  cenário da diferença (grade, dispersão 0,001, extirpação 0,02, alfa contando sítios vazios):
  1,51 com Poisson, 1,61 com multinomial e 1,61 no mcomsimr. O multinomial custa o mesmo em
  dispersão baixa e até 20 vezes mais em dispersão alta (10 x 10 sítios, 50 espécies).
- **Erros no artigo de Thompson et al. (2020).** Tabela 1 (beta), faixa de dispersão,
  conversão do ambiente em postos e metade dos imigrantes no cenário de
  competição-colonização não batem com o código. Vale citar isso no artigo do pacote.
- **Dependências do hub.** Decidido em 2026-11-07: `Imports` só com `ape` e pacotes base
  (`graphics`, `grDevices`, `parallel`, `stats`, `utils`); todo o resto em `Suggests`, chamado
  com `requireNamespace()`. O `purrr`, que só dava a barra de progresso, foi para `Suggests`.
- **Código em MATLAB.** Os suplementos não trazem licença. Decidido em 2026-11-21 (fase 2 do
  rOpenSci): pedir aos autores (Peres-Neto e Leibold) permissão para distribuir as traduções
  sob MIT e perguntar se querem constar em `Authors@R` (`ctb` ou `cph`), e corrigir já a
  redação (`DESCRIPTION` dizia que tudo vinha das equações publicadas). E-mail enviado em
  2026-11-21; aguardando resposta. Se recusarem ou não
  responderem, reescrever `node_test()`, `node_contrasts()`, `phylo_vector()` e
  `trait_partition()` a partir dos artigos e apêndices, com o Octave só para conferir.
  `food_web_dynamics()` foi escrito a partir das equações; o MATLAB de McPeek só foi rodado.
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
- **Validar `simulate_logistic()` com PGLMM.** Feito em 2026-11-07, em
  `data-raw/validation_pglmm.R` (phyr 1.1.3, só no script; árvore balanceada de 32 espécies,
  40 sítios, 100 réplicas por cenário, 900 ajustes sem falha). Modelo
  `y ~ x + (1 | sp) + (x | sp__)`, binomial. Sem a redução por sítio, o PGLMM recupera a
  simulação: fração filogenética da variância das inclinações 0,00, 0,48 e 1,00 para lambda
  0, 0,5 e 1, igual à obtida das inclinações verdadeiras; teste da inclinação filogenética com
  5-6% de rejeição em lambda 0, 71% em 0,5 e 100% em 1; variâncias 10-15% abaixo, viés
  conhecido do ajuste binário. Com a redução padrão U(0,5; 1), que o modelo não tem, as
  variâncias caem à metade, mas a fração e as rejeições não mudam; um intercepto aleatório de
  sítio não corrige. `slope_sd` é o desvio amostral das inclinações (padronizadas por
  `scale()`), não o do processo: com lambda 1 nesta árvore, a variância do processo é 1,48.
  Isso ficou explícito na ajuda.
- **treeshapesim.** O artigo de integração com o comdynamite usa `species_traits()`,
  `competition_matrix()`, `landscape_implicit()` e `simulate_metacommunity()`: mudar
  assinatura ou saída delas quebra o site do treeshapesim. Quando `extinct_tips()` (antes
  `prune_tips()`, renomeada no PR #8 de lá) estiver na `main` de lá, pôr o link do artigo no README ou em `evaluating-methods`. O
  comdynamite fica só em `Config/Needs/website` no treeshapesim, sem dependência circular.
- **Sessão anterior.** A sessão que começou o pacote não foi avisada por aqui de que o
  desenvolvimento mudou de lugar. Commits novos dela vão divergir do `main`.
- **pcps-sim.** Se a submissão for duplo-cega, avaliar se o `NEWS.md` pode manter "Code
  from the pcps-sim project". A vignette `evaluating-methods` compara permutações
  livres e do toro, tema próximo do piloto de validade sob H0 daquele projeto.
- **Resultados pré-calculados.** Os `.rds` em `vignettes/` vêm de
  `data-raw/vignette_ecology.R` e `data-raw/vignette_method.R`. Se o modelo mudar,
  é preciso rodar os scripts de novo e revisar os números citados no texto.
- **Páginas órfãs do site.** Resolvido em 2026-11-07: o deploy passou a `clean: true` (com
  `clean-exclude: .nojekyll`), que apaga do `gh-pages` o que o site não gera mais. Na
  comparação com um site construído localmente, sobravam só 5 figuras antigas de artigos.
- **`data-raw/hex_dynamite.R`.** Removido em 2026-11-07; o logo atual (opção F1) vem de
  `data-raw/hex_sticker.R`.

## Referências citadas neste roadmap

- Desdevises, Y., Legendre, P., Azouzi, L., & Morand, S. (2003). *Evolution*, 57, 2647–2652. https://doi.org/10.1111/j.0014-3820.2003.tb01508.x
- Desjardins-Proulx, P., & Gravel, D. (2012a). *The American Naturalist*, 179, 137–144. https://doi.org/10.1086/663196
- Desjardins-Proulx, P., & Gravel, D. (2012b). *Ecology and Evolution*, 2, 1781–1790. https://doi.org/10.1002/ece3.292
- Drury, J., Clavel, J., Manceau, M., & Morlon, H. (2016). *Systematic Biology*, 65, 700–710. https://doi.org/10.1093/sysbio/syw020
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
- Hubert, N., Calcagno, V., Etienne, R. S., & Mouquet, N. (2015). *Ecology Letters*, 18, 864–881. https://doi.org/10.1111/ele.12458
- Harvey, M. G., & Rabosky, D. L. (2017). *Methods in Ecology and Evolution*, 9, 984–993. https://doi.org/10.1111/2041-210x.12949
- Ives, A. R., & Helmus, M. R. (2011). *Ecological Monographs*, 81, 511–525. https://doi.org/10.1890/10-1264.1
- Legendre, P. (2008). *Journal of Plant Ecology*, 1, 3–8. https://doi.org/10.1093/jpe/rtm001
- Januario, M., Pinsky, M. L., & Rabosky, D. L. (2025). *Ecology Letters*, 28, e70021. https://doi.org/10.1111/ele.70021
- Leibold, M. A., Govaert, L., Loeuille, N., De Meester, L., & Urban, M. C. (2022). *Annual Review of Ecology, Evolution, and Systematics*, 53, 299–326. https://doi.org/10.1146/annurev-ecolsys-102220-024934
- Leibold, M. A., Urban, M. C., De Meester, L., Klausmeier, C. A., & Vanoverbeke, J. (2019). *Proceedings of the National Academy of Sciences*, 116, 2612–2617. https://doi.org/10.1073/pnas.1808615116
- Leibold, M. A., et al. (2010). *Ecology Letters*, 13, 1290–1299. https://doi.org/10.1111/j.1461-0248.2010.01523.x
- Matzke, N. J. (2014). *Systematic Biology*, 63, 951–970. https://doi.org/10.1093/sysbio/syu056
- McPeek, M. A. (2017). *Evolutionary Community Ecology*. Monographs in Population Biology. Princeton University Press. https://doi.org/10.23943/princeton/9780691088778.001.0001
- Morlon, H., et al. (2016). *Methods in Ecology and Evolution*, 7, 589–597. https://doi.org/10.1111/2041-210x.12526
- O'Meara, B. C., Ané, C., Sanderson, M. J., & Wainwright, P. C. (2006). *Evolution*, 60, 922–933. https://doi.org/10.1111/j.0014-3820.2006.tb01171.x
- Nakamura, G., et al. (2023). *Journal of Biogeography*, 51, 533–543. https://doi.org/10.1111/jbi.14610
- Pontarp, M., et al. (2019). *Methods in Ecology and Evolution*, 10, 450–460. https://doi.org/10.1111/2041-210x.13129
- Prangle, D., et al. (2014). *Statistical Applications in Genetics and Molecular Biology*, 13, 67–82. https://doi.org/10.1515/sagmb-2013-0012
- Pudlo, P., et al. (2015). *Bioinformatics*, 32, 859–866. https://doi.org/10.1093/bioinformatics/btv684
- Raynal, L., et al. (2019). *Bioinformatics*, 35, 1720–1728. https://doi.org/10.1093/bioinformatics/bty867 (verificada no Crossref em 2026-11-11)
- Roberts, D. R., et al. (2017). *Ecography*, 40, 913–929. https://doi.org/10.1111/ecog.02881 (verificada no Crossref em 2026-11-11)
- Rabosky, D. L., & Goldberg, E. E. (2015). *Systematic Biology*, 64, 340–355. https://doi.org/10.1093/sysbio/syu131
- Rabosky, D. L., & Goldberg, E. E. (2017). *Evolution*, 71, 1432–1442. https://doi.org/10.1111/evo.13227
- McArdle, B. H., & Anderson, M. J. (2001). *Ecology*, 82, 290–297. https://doi.org/10.1890/0012-9658(2001)082[0290:FMMTCD]2.0.CO;2
- Peres-Neto, P. R., & Kembel, S. W. (2015). *Plant Ecology*, 216, 709–724. https://doi.org/10.1007/s11258-014-0405-0
- Peres-Neto, P. R., et al. (2012). *Ecology*, 93, S14–S30. https://doi.org/10.1890/11-0494.1
- Pigot, A. L., & Etienne, R. S. (2015). *Ecology Letters*, 18, 153–163. https://doi.org/10.1111/ele.12395
- Robert, C. P., et al. (2011). *PNAS*, 108, 15112–15117. https://doi.org/10.1073/pnas.1102900108
- Ruffley, M., et al. (2019). *Ecology and Evolution*, 9, 13218–13230. https://doi.org/10.1002/ece3.5773
- van der Plas, F., et al. (2015). *Ecology*, 96, 1502–1515. https://doi.org/10.1890/14-0454.1
- Urban, M. C., Strauss, S. Y., Pelletier, F., Palkovacs, E. P., Leibold, M. A., Hendry, A. P., De Meester, L., Carlson, S. M., Angert, A. L., & Giery, S. T. (2020). *Proceedings of the National Academy of Sciences*, 117, 17482–17490. https://doi.org/10.1073/pnas.1918960117
- Violle, C., Enquist, B. J., McGill, B. J., Jiang, L., Albert, C. H., Hulshof, C., Jung, V., & Messier, J. (2012). *Trends in Ecology & Evolution*, 27, 244–252. https://doi.org/10.1016/j.tree.2011.11.014
- Vehtari, A., et al. (2017). *Statistics and Computing*, 27, 1413–1432. https://doi.org/10.1007/s11222-016-9696-4
- Wood, S. N. (2010). *Nature*, 466, 1102–1104. https://doi.org/10.1038/nature09319

Conferidas no Crossref em 2026-11-10 (item do `pcps.curve()`).

- Debastiani, V. J., & Duarte, L. D. S. (2014). *Frontiers of Biogeography*, 6, 144–148. https://doi.org/10.21425/F5FBG22943
- Diniz-Filho, J. A. F., Rangel, T. F., Santos, T., & Bini, L. M. (2012). *Evolution*, 66, 1079–1090. https://doi.org/10.1111/j.1558-5646.2011.01499.x
- Duarte, L. D. S., Debastiani, V. J., Carlucci, M. B., & Diniz-Filho, J. A. F. (2018). *Ecology*, 99, 385–398. https://doi.org/10.1002/ecy.2081
- Peres-Neto, P. R., Dray, S., & ter Braak, C. J. F. (2017). *Ecography*, 40, 806–816. https://doi.org/10.1111/ecog.02302
- Zelený, D. (2018). *Journal of Vegetation Science*, 29, 953–966. https://doi.org/10.1111/jvs.12688
