# Plano: especiação acoplada à metacomunidade

Objetivo: fechar o limite registrado em `dev/notas-artigo-taxas.md` (ponte Leibold–Januario).
Hoje o comdynamite simula a metacomunidade sobre um pool fixo de espécies e uma árvore dada;
a especiação não depende da dinâmica. Com populações locais que divergem e viram espécies, a
árvore passa a emergir da metacomunidade, e dá para perguntar com verdade conhecida se a
persistência metapopulacional prediz a especiação (Januario et al. 2025) e se a montagem
deixa rastro nas taxas de evolução de traits entre espécies (Leibold et al. 2019, 2022).

Referências conferidas no Crossref ao fim deste arquivo. O que não é afirmação de uma fonte é
proposta nossa.

## Situação atual do código

- `mc_dynamics()` (em `R/dynamics.R`) integra abundâncias de k espécies em m sítios com
  crescimento dependente do ambiente (`r`), competição (`alpha`), dispersão, imigração do
  pool e extirpação. k é fixo: as matrizes têm uma linha por espécie desde o início.
- Os traits são por espécie, não por população; não há evolução durante a dinâmica.
- `trait_mk()`, `history_rates()` e `trait_rates()` simulam traits e taxas sobre uma árvore
  dada. `food_web_dynamics()` tem evolução de traits médios por gradiente, mas sem espaço.

## Modelo proposto

Segue a família de modelos mecanísticos eco-evolutivos em que especiação, dispersão e
seleção ambiental geram a diversidade regional, revisada por Pontarp et al. (2019) e
implementada de forma geral no gen3sis (Hagen et al. 2021). O que o comdynamite acrescenta é a
ligação com os métodos de metacomunidade e de taxas por ramo do pacote, com a verdade
conhecida em cada etapa.

### 1. Populações com trait local

- Unidade: população = espécie × sítio, com abundância (já existe) e trait médio local z_is
  (novo).
- Seleção: z_is anda pelo gradiente de ln W no sítio, como em `food_web_dynamics()`
  (dz/dt = V ∂ln W/∂z, McPeek 2017, eq. 3.7), com W dado pelo ajuste ao ambiente local e,
  opcionalmente, pela competição dependente do trait.
- Fluxo gênico: depois da dispersão, z_is vira a média ponderada pelas abundâncias de
  residentes e imigrantes. Com dispersão alta isso afoga a adaptação local (gene swamping,
  Leibold et al. 2022); com dispersão baixa permite os efeitos de prioridade mediados pela
  evolução de Leibold et al. (2019).
- Variância genética V fixa por espécie na primeira versão; genética explícita (locos, como no
  modelo baseado em indivíduos de Leibold et al. 2019) fica como extensão.

### 2. Divergência e especiação

Modos, escolhidos por argumento, porque a regra de especiação condiciona os resultados e
precisa ser variada nos experimentos:

- `"allopatric"` (modo principal): uma matriz de divergência entre populações coespecíficas
  cresce a cada passo sem fluxo gênico entre elas e é zerada (ou reduzida) quando há fluxo.
  Quando a divergência entre grupos de populações passa de um limiar, o grupo vira espécie
  nova. É a abordagem do gen3sis (Hagen et al. 2021), reimplementada a partir da descrição
  publicada.
- `"ecological"`: especiação quando a diferença de trait entre grupos de populações passa de
  um limiar, com ou sem exigência de isolamento.
- `"protracted"`: espécies incipientes que só viram espécies boas depois de um tempo τ, como
  em Rosindell et al. (2010) e Etienne & Rosindell (2012); permite comparar com o modelo PBD.
- `"point"`: mutação pontual neutra (cada indivíduo novo vira espécie nova com probabilidade
  ν), para o limite neutro de teste.

### 3. Registro filogenético e da história

- Cada especiação registra tempo, espécie mãe, espécies filhas e as populações de cada uma.
  No fim, montar um `phylo` completo (com extintas) e o reconstruído (só sobreviventes, via
  `ape::drop.fossil()` ou equivalente).
- Extinção global quando a abundância regional de uma espécie zera.
- A distribuição de cada linhagem ao longo do tempo (sítios ou tipos de hábitat ocupados) vira
  um `simmap`, pronto para `history_rates()` e `as_herodotools()`.
- Séries por espécie: ocupação, abundância regional e as duas métricas de persistência de
  Januario et al. (2025), ocupação de equilíbrio e resiliência (conferir a definição exata no
  artigo antes de implementar).
- Traits por população e médias por espécie, para os estimadores de taxa por ramo e para as
  estatísticas T (`cati::Tstats()`).

### 4. Interface (decidir na rodada de UX)

- Opção A: função nova, por exemplo `simulate_diversification(landscape, ...)`, que reaproveita
  as partes internas de `mc_dynamics()` (dispersão, imigração, competição) e devolve uma
  classe própria com `tree`, `tree_full`, `history` (simmap), `abundance`, `pop_traits`,
  `events` e `persistence`.
- Opção B: argumento `speciation = ` em `simulate_metacommunity()`. Mais simples para o
  usuário, mas complica uma função que hoje assume k fixo.
- Preferência inicial: A, com um `as_metacommunity()` que converta a saída para os métodos que
  já existem.

### 5. Implementação

- R base, sem dependência nova. O número de espécies cresce: guardar abundâncias e traits em
  matrizes que crescem por blocos (dobrar a capacidade) e um teto `max_species` com erro
  claro.
- Medir o custo cedo (como em `simulation_cost()`), porque a matriz de divergência é
  quadrática no número de populações por espécie. Se ficar lento demais, discutir Rcpp antes
  de adotar, porque é dependência nova.

## Validação

1. Limite neutro (`"point"`, sem seleção): o número de espécies e as abundâncias devem bater
   com as expectativas da teoria neutra para o θ correspondente. Antes, levantar quais
   expressões analíticas usar (Rosindell et al. 2010 dá resultados para o caso prolongado).
2. Modo prolongado contra o PBD: as árvores reconstruídas devem ter a desaceleração aparente
   descrita por Etienne & Rosindell (2012). Ajustar o PBD (pacote `PBD`, conferir a licença) em
   `data-raw` e ver se recupera as taxas.
3. gen3sis (CRAN, GPL; só em `data-raw`): mesma paisagem e parâmetros equivalentes no modo
   alopátrico; comparar riqueza, tamanho de clado e forma da árvore.
4. Leibold et al. (2019): sem especiação e com duas espécies, reproduzir a neutralidade
   regional qualitativamente (convergência das distribuições de trait e das ocupações) nos
   eixos de extinção de manchas e dispersão da Fig. 5.
5. Consistência interna: árvore ultramétrica no reconstruído, número de eventos igual ao de
   nós, simmap coerente com as ocupações, testes de que sem especiação o resultado é o de
   `mc_dynamics()` com a mesma semente.

## Experimentos que isso abre

- Teste de Januario com verdade conhecida: persistência por espécie contra a taxa de
  especiação nas pontas da árvore reconstruída (por exemplo a estatística DR de Jetz et al.
  2012), variando dispersão e extinção de manchas. Pergunta: o desacoplamento aparece mesmo
  quando a persistência causa a especiação no modelo?
- Neutralidade regional e taxas de trait: taxas por ramo verdadeiras e estimadas com médias de
  espécies, contra a divergência entre populações e as estatísticas T.
- Ligação com os cenários 1, 2, 5 e 8 das notas, agora com a árvore gerada pela própria
  metacomunidade.

## Fases

1. Traits por população com seleção e fluxo gênico em `mc_dynamics()`, sem especiação;
   validação 4.
2. Divergência, os quatro modos de especiação e o registro filogenético; validações 1, 2 e 5.
3. Saídas (simmap, persistência, traits por população) e integração com `history_rates()`,
   `trait_rates()` e as estatísticas T.
4. Comparação com o gen3sis (validação 3) e medida de custo.
5. Experimentos e o artigo.

## Riscos e decisões em aberto

- A regra de especiação decide muito do resultado; por isso os vários modos e a obrigação de
  reportar a sensibilidade.
- Escalas: a dinâmica roda em passos ecológicos e a especiação precisa de muitos passos.
  Januario et al. (2025) atribuem o desacoplamento justamente à diferença entre décadas e
  milhões de anos; o modelo comprime essas escalas e o artigo tem de dizer isso.
- Custo computacional com muitas espécies e sítios.
- Decisões para o Diogo: nome e forma da interface (A ou B); quais modos entram na primeira
  versão; número da versão (proposta: bloco próprio, número a definir, depois da v0.7.0).

## Referências

- Etienne, R. S., & Rosindell, J. (2012). *Systematic Biology*, 61, 204 (página final não conferida).
  https://doi.org/10.1093/sysbio/syr091
- Hagen, O., Flück, B., Fopp, F., Cabral, J. S., Hartig, F., Pontarp, M., Rangel, T. F., &
  Pellissier, L. (2021). *PLOS Biology*, 19, e3001340.
  https://doi.org/10.1371/journal.pbio.3001340
- Januario, M., Pinsky, M. L., & Rabosky, D. L. (2025). *Ecology Letters*, 28, e70021.
  https://doi.org/10.1111/ele.70021
- Jetz, W., Thomas, G. H., Joy, J. B., Hartmann, K., & Mooers, A. O. (2012). *Nature*, 491,
  444–448. https://doi.org/10.1038/nature11631
- Leibold, M. A., Urban, M. C., De Meester, L., Klausmeier, C. A., & Vanoverbeke, J. (2019).
  *PNAS*, 116, 2612–2617. https://doi.org/10.1073/pnas.1808615116
- Leibold, M. A., Govaert, L., Loeuille, N., De Meester, L., & Urban, M. C. (2022). *Annual
  Review of Ecology, Evolution, and Systematics*, 53, 299–326.
  https://doi.org/10.1146/annurev-ecolsys-102220-024934
- McPeek, M. A. (2017). *Evolutionary Community Ecology*. Princeton University Press.
  https://doi.org/10.23943/princeton/9780691088778.001.0001
- Pontarp, M., et al. (2019). *Trends in Ecology & Evolution*, 34, 211–223.
  https://doi.org/10.1016/j.tree.2018.11.009
- Rosindell, J., Cornell, S. J., Hubbell, S. P., & Etienne, R. S. (2010). *Ecology Letters*,
  13, 716–727. https://doi.org/10.1111/j.1461-0248.2010.01463.x
