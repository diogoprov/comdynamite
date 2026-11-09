# Plano: artigo sobre montagem local e regimes de diversificação (DAMOCLES)

Escrito em 2026-11-09. Documento interno de trabalho; nada aqui foi rodado ainda. Cada
bloco de simulação e de análise empírica só começa depois de aprovado. Não confundir com
`dev/notas-artigo-taxas.md`, que trata de taxas de evolução de traits.

## 1. Pergunta e enquadramento

**Pergunta.** As taxas de colonização (γ) e de extinção local (μ) das comunidades locais
diferem entre regimes de diversificação, além do que a própria história de especiação já
explica?

**Por que agora.** Reijenga et al. (2026) estenderam o DAMOCLES para vários clados ao mesmo
tempo. Eles mostraram que clados de diversificação rápida têm proporcionalmente menos
espécies em simpatria e atribuíram isso só à história de especiação: γ e μ eram constantes
entre clados, e o modelo global venceu o modelo por família por AIC. Ninguém testou se γ e μ
mudam *dentro* de uma árvore, acompanhando regimes ou taxas de diversificação estimados por
outro método. É essa a ponte que o pacote passou a fazer (PRs #39, #40, #44 e #56), com:

- regimes e covariáveis de ramo em `damocles_fit()`;
- ajustes sobre configurações do posterior em `damocles_samples()`;
- teste contra partições aleatórias em `damocles_partition_test()`;
- taxas de ponta e de ramo com `dr_rates()` e os exemplos de BAMM e ClaDS.

A busca bibliográfica de 2026-11-09 (OpenAlex, Scholar Gateway e Undermind) não encontrou
precedente para essa integração.

**Precedente por passo:**

- o ajuste compartilhado já tem precedente direto (Reijenga et al. 2026);
- a simpatria tem análogos em DEC e GeoSSE;
- LRT e intervalos por perfil são inferência padrão;
- a novidade está em regimes e covariáveis vindos de BAMM e ClaDS, na propagação do
  posterior, no teste de partições e na diluição por regressão.

**Formato proposto (a decidir).**

- *Artigo 1, de métodos (MEE ou Syst Biol):* o modelo, as simulações da seção 3 e
  Benício et al. (2021) como exemplo trabalhado.
- *Artigo 2, biológico (AmNat; Ecol Lett como alternativa):* a pergunta na escala da Mata
  Atlântica, com o ATLANTIC AMPHIBIANS e vários clados.

Opção alternativa: um artigo único para a AmNat. O risco é ter revisores de métodos e de
ecologia cobrando coisas diferentes.

## 2. Dados

Tudo em `data-raw/` no checkout principal, sem commit. A pasta `data-raw` está no
`.Rbuildignore`. Licenças a respeitar antes de qualquer redistribuição.

| Conjunto | Onde | Conteúdo | Situação |
|---|---|---|---|
| Benício et al. (2021) | `data-raw/benicio2021/` (Dryad, doi:10.5061/dryad.sf7m0cg65) | presença e ausência de anuros em 14 áreas protegidas (7 no interior, 7 na costa), coordenadas, `geo_frogs.txt`, modos reprodutivos, scripts (GeoSSE, BioGeoBEARS) | baixado |
| Árvores de Benício | mesma pasta | `tree_pyron_112_species` e `tree_feng_112_species` (casam 112/112 com a matriz); `tree_pyron_248_species` e `tree_feng_248_species` (NEXUS, 97 árvores), com 95–96/112 por sinônimos (*Hypsiboas* → *Boana*, *Scinax*) | baixado |
| Portik et al. (2023) | `data-raw/portik2023/` (OSF 3cbsh) | árvore datada por TreePL, 5.326 pontas, binária e ultramétrica, mais 100 bootstraps datados | baixado; o OSF só diz "copyright Daniel Portik, 2023", sem licença aberta: **não redistribuir** |
| Vasquez-Restrepo (2024, preprint) | Figshare doi:10.6084/m9.figshare.24894363 (CC BY 4.0) | BAMM em 100 árvores totalmente amostradas de anfíbios; DR; presença e ausência da IUCN | não baixado (`1. Diversification.zip`, 6,7 GB) |
| ATLANTIC AMPHIBIANS (Vancine et al. 2018) | `af-frogs-climate-rarity/data/raw/atlantic_amphibians/`; também no `atlanticr` | 1.163 sítios, 17.626 registros | disponível; esforço amostral heterogêneo |

### 2.1 Casamento de nomes com a árvore do Portik (feito em 2026-11-09)

- Normalizando aspas e `_` final e trocando *Hypsiboas* por *Boana*:
  - 91/112 espécies da matriz estão na árvore do Portik;
  - 214/248 do pool amplo estão na árvore do Portik.
- As faltas são de quatro tipos:
  1. **concordância de gênero em *Boana*:** *albomarginata*, *bischoffi*, *pulchella*,
     *semilineata* etc. Resolve por tabela;
  2. **erros de grafia:** *Bokermanohyla*, *Aparasphenedon*, *Lysapus*, *Crossodactylus
     schmdti*. Resolve por tabela;
  3. **espécies não descritas ou "aff.":** *Adenomera* sp. J, *Aplastodiscus* sp. 5,
     *Fritziana* aff. *fissilis*, *Megaelosia* aff. *boticariana*, *Scinax* aff. *hayii*,
     *Dendrophryniscus* sp. Precisam ser enxertadas junto à congênere mais próxima, com a
     árvore de Benício como guia;
  4. **linhagens de *Ischnocnema guentheri* (C1, C3, C4) e *I. henseli*:** o Portik não
     separa essas linhagens. Decidir entre colapsar em *I. guentheri* e enxertar.
- Entregável: `data-raw/benicio2021_taxonomy.csv` (nome original → nome na árvore → ação)
  e um script que constrói a árvore do pool a partir do Portik e dos 100 bootstraps.

### 2.2 Definição do pool regional (decisão central)

- Com a árvore de 112 espécies, o pool é a união das 14 comunidades: nenhuma espécie está
  ausente de todos os sítios. Isso condiciona o DAMOCLES a linhagens que colonizaram algum
  lugar e tende a superestimar γ.
- As árvores de 248 espécies sugerem um pool mais amplo, com 136 espécies ausentes das 14
  áreas. **Falta confirmar o que elas representam** (README do Dryad, ou perguntar ao
  Fernando ou ao Ronildo).
- Plano:
  - pool de 248, harmonizado e podado do Portik, como análise principal;
  - pool de 112 como sensibilidade;
  - a simulação S6 quantifica o viés esperado.
- Costa e interior: o DAMOCLES assume um único pool. Duas formas, a comparar:
  - (a) pool combinado com região como regime de γ e μ nos ramos terminais, o que exige
    estender `regimes` para depender da comunidade, não só do ramo (ver 3.8);
  - (b) ajustes separados por região, comparados por LRT contra o modelo compartilhado.

## 3. Lacunas metodológicas e simulações

Todas em `data-raw/`, um script por bloco, com:

- `RNGkind("L'Ecuyer-CMRG")`, `CORES` por variável de ambiente e `.rds` por cenário
  (pontos de retomada);
- `message()` de progresso por cenário;
- nenhum `repeat` sem limite de tentativas (lição da seção 7 da validação, que travou por
  cerca de 6 h).

Árvores simuladas com `treeshapesim` ou com processos de nascimento e morte com mudanças de
taxa (pacotes GPL só em `data-raw/`).

Os tempos abaixo são **estimativas grosseiras**, a medir num piloto antes de cada bloco.

### S1. Erro tipo I com diversificação heterogênea (prioridade máxima)

- **Pergunta:** o teste de regime rejeita γ constante só porque a árvore tem mudanças de
  taxa de especiação?
- **Desenho:**
  - árvores com 1–3 mudanças de λ (razão 2–4×), com 112 e 250 pontas;
  - comunidades simuladas com γ e μ constantes;
  - n = 14 e n = 40 comunidades.
- **Regimes testados:**
  - (i) os verdadeiros;
  - (ii) os estimados, com DR por ramo (barato), ClaDS (moderado) e BAMM num subconjunto
    de cerca de 50 réplicas (caro).
- **Saída:** proporção de rejeição pelo χ² e pelas partições aleatórias, com IC binomial.
  Meta: 5%.
- Estimativa: 1–2 dias de CPU com 6 núcleos, sem BAMM.

### S2. Poder

- Grade:
  - n ∈ {5, 10, 14, 20, 40};
  - pontas ∈ {50, 112, 250};
  - razão de γ entre regimes ∈ {1,5; 2; 3};
  - efeito em γ, em μ ou em ambos;
  - fração da árvore no regime focal (10%, 25%, 50%).
- Inclui o caso realista de Benício: n = 14 e pool de 112 ou 248.
- Saída: curvas de poder. É ela que diz se n = 14 sustenta o caso focal sozinho.

### S3. Integração do posterior (o ponto metodologicamente aberto)

Hoje `damocles_samples()` resume por "proporção de amostras com p < 0,05" e β mediano. A
diluição já mostrada com o ClaDS (β mediano das amostras 0,22; MAP 1,64; verdadeiro 2)
torna esse resumo inadequado para covariáveis contínuas.

Estimadores a comparar:

- (a) mediana por amostra, o atual;
- (b) taxas médias do posterior por ramo como covariável;
- (c) **verossimilhança integrada sobre configurações**: L̄(θ) = (1/K) Σ_k L(θ | config_k),
  maximizada uma vez. É natural para regimes discretos do BAMM e propaga a incerteza das
  configurações em vez de resumi-la. A ressalva a checar é que as configurações do BAMM
  foram amostradas condicionalmente à árvore, e não aos dados de comunidade;
- (d) correção de erro nas variáveis (SIMEX ou calibração de regressão) para covariáveis
  contínuas, usando a variância posterior por ramo.

Métricas: viés e RMSE de β, cobertura dos intervalos por perfil, erro tipo I e poder.

Entregável provável: um argumento `integrate = TRUE` (ou uma função nova) que implemente
(c), e uma recomendação documentada para (b) ou (d). Isso vira código no pacote, com testes,
num PR próprio.

### S4. Calibração do teste de partições aleatórias

- Indício atual: 0,10 de rejeição no cenário nulo com 40 réplicas e 0,17 num teste menor.
- 500 réplicas do nulo (com e sem diversificação heterogênea) e variação de `tolerance`,
  `max_overlap` e `n_rand`.
- Se o teste for liberal, avaliar:
  - mid-p;
  - restringir o sorteio a clados com idade de coroa parecida, além do tamanho;
  - mais partições por réplica.

### S5. Dependência espacial entre comunidades

- O ajuste compartilhado trata as comunidades como independentes, dada a árvore. Sítios
  próximos compartilham história de colonização, e o LRT tende a ficar anticonservador.
- **Simulação:** comunidades em blocos que compartilham os estados até um tempo de
  separação t_s, e só então evoluem independentemente. t_s = 0 equivale a comunidades
  independentes.
- **Avaliar:**
  - inflação do LRT em função de t_s;
  - bootstrap em blocos de sítios (por região, ou por agrupamento espacial) como correção;
  - partições aleatórias com reamostragem em blocos.
- Com Benício, os blocos naturais são as duas regiões.

### S6. Pool mal especificado

- Simular com um pool verdadeiro de 250 espécies e ajustar com:
  - (i) o pool verdadeiro;
  - (ii) a união das espécies observadas (como a árvore de 112).
- Medir o viés em γ, μ e β, e o efeito no teste de regime.

### S7. Incerteza de topologia e datação

- Caso empírico: repetir o ajuste principal nos 100 bootstraps datados do Portik e nas
  97 árvores da Feng.
- Simulação: só um teste rápido de que a variação entre árvores não infla o erro tipo I.

### 3.8 Extensões de código que o plano exige

1. Verossimilhança integrada sobre configurações (S3c).
2. Bootstrap em blocos de comunidades (S5), em `damocles_partition_test()` ou numa função
   de reamostragem.
3. Regime que depende da comunidade (costa e interior) nos ramos terminais (2.2a), se a
   opção (a) for escolhida.
4. Funções auxiliares de casamento de nomes e enxerto ficam em `data-raw/`, não no pacote.

## 4. Caso de estudo: ordem das análises

- **E1. Taxonomia e árvore do pool:** a tabela de 2.1, a poda do Portik até 248 (e até
  112), o enxerto das espécies faltantes e os 100 bootstraps podados.
- **E2. Regimes e taxas:**
  - BAMM do Vasquez-Restrepo, se o zip de 6,7 GB tiver os dados de eventos. Confirmar
    listando o zip sem baixar tudo, e baixar só o que for de anfíbios;
  - ClaDS próprio sobre o Portik podado (RPANDA, em `data-raw/`);
  - DR com `dr_rates()`;
  - GeoSSE de Benício et al. (2021) como referência da diferença costa contra interior.
- **E3. Ajuste base:** `damocles_fit(shared = TRUE)` com γ e μ constantes, pools de 248 e
  112, e `sympatry` livre e fixa em 0.
- **E4. Regimes e covariáveis:**
  - LRT e intervalos por perfil;
  - `damocles_samples()` sobre as configurações do posterior, com o estimador vencedor
    de S3;
  - `damocles_partition_test()` com a melhor configuração.
- **E5. Robustez:**
  - bootstrap em blocos (S5);
  - 100 árvores do Portik e 97 da Feng (S7);
  - pool de 112 contra 248;
  - com e sem as espécies enxertadas.
- **E6. Escala (artigo 2):** o ATLANTIC AMPHIBIANS filtrado por esforço e por qualidade
  da coordenada, com ajustes por clado (famílias ou gêneros com pelo menos cerca de 20
  espécies) e regimes dentro de cada clado. Comparar com os modelos "global" e "por clado"
  de Reijenga et al. (2026).

## 5. Cronograma (estimativa)

| Etapa | Duração |
|---|---|
| E1 (taxonomia e árvore) + piloto de tempo de computação | 1 semana |
| S1, S2, S4, S6 (em paralelo, em background) | 3–5 semanas |
| S3 + código da verossimilhança integrada + S5 e bootstrap em blocos | 2–3 semanas |
| E2–E5 (caso de estudo) | 2–3 semanas |
| Redação do artigo 1 | 4–6 semanas |
| E6 (artigo 2) | a planejar depois do artigo 1 |

## 6. Decisões pendentes (do Diogo)

1. Artigo 1 de métodos e artigo 2 biológico, ou artigo único para a AmNat.
2. O que é o pool de 248 espécies de Benício.
3. Linhagens de *Ischnocnema guentheri*: colapsar ou enxertar.
4. Se baixamos as partes de anfíbios do Figshare do Vasquez-Restrepo (até 6,7 GB).
5. Autoria e contato com Etienne e Pigot (Reijenga et al. 2026 é do mesmo grupo).

## 7. Riscos

- **Poder com n = 14 pode ser baixo** (S2). Nesse caso, Benício vira ilustração e a
  inferência pesa no ATLANTIC AMPHIBIANS.
- **Super-árvores com espécies imputadas** (Vasquez-Restrepo) distorcem taxas perto das
  pontas. Mitigação: comparar com o Portik, que é molecular.
- **Confusão entre barreira biogeográfica e dinâmica local** quando costa e interior
  compartilham o pool (2.2).
- **Esforço heterogêneo no ATLANTIC AMPHIBIANS:** ausências falsas inflam μ. Pode exigir
  um modelo de detecção, o que está fora do escopo do artigo 1.

## Referências (verificadas no Crossref, OpenAlex ou servidor de preprints em 2026-11-09)

- Benício, R., Provete, D. B., Lyra, M., Heino, J., Haddad, C. F. B., Rossa-Feres, D. C., &
  da Silva, F. R. (2021). Differential speciation rates, colonization time and niche
  conservatism affect community assembly across adjacent biogeographical regions.
  *Journal of Biogeography*, 48, 2211–2225. doi:10.1111/jbi.14145
- Pigot, A. L., & Etienne, R. S. (2015). A new dynamic null model for phylogenetic
  community structure. *Ecology Letters*, 18, 153–163. doi:10.1111/ele.12395
- Portik, D. M., Streicher, J. W., Blackburn, D. C., Moen, D. S., et al. (2023). Redefining
  possible: combining phylogenomic and supersparse data in frogs. *Molecular Biology and
  Evolution*, 40, msad109. doi:10.1093/molbev/msad109
- Reijenga, B. R., Etienne, R. S., Murrell, D. J., & Pigot, A. L. (2026). Speciation
  history shapes patterns of assemblage species richness in birds. *Ecology Letters*, 29,
  e70432. doi:10.1111/ele.70432
- Vancine, M. H., Duarte, K. S., de Souza, Y. S., et al. (2018). ATLANTIC AMPHIBIANS: a
  data set of amphibian communities from the Atlantic Forests of South America. *Ecology*,
  99, 1692. doi:10.1002/ecy.2392
- Vasquez-Restrepo, J. D. (2024). A dataset containing speciation rates, uncertainty, and
  presence-absence matrices for more than 34,000 vertebrate species. bioRxiv (preprint,
  sem revisão por pares). doi:10.1101/2024.04.09.588748
- Para o restante do contexto metodológico (STRAPP, ES-sim, Title & Rabosky 2019, ClaDS,
  BAMM), ver a documentação de `damocles_samples()` e `dr_rates()`.
