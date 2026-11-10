# Notas para o artigo sobre taxas por ramo

Leituras para o artigo (vignette) sobre taxas de evolução de traits que variam entre ramos e o
papel da montagem de metacomunidades. Distingue o que cada fonte afirma da nossa interpretação.
Citações literais conferidas contra o PDF ou o texto extraído. Referências completas na seção
"Referências" do `ROADMAP.md`.

## Pergunta do artigo

A montagem da metacomunidade molda a evolução dos traits e deixa rastro nas taxas estimadas por
ramo? `trait_rates()` e `trait_mk()` dão a verdade conhecida; a validação em
`data-raw/validation_rates.R` e `data-raw/validation_bamm.R` mostra quanto os estimadores
recuperam. Limite central: com uma realização por ramo, a taxa de cada ramo é pouco
identificável (correlação entre mudança verdadeira e taxa verdadeira de 0,29 a 0,40).

## Hubert et al. (2015)

- Revisão dos modelos de especiação em metacomunidades (neutros, de nicho, de trade-offs de
  história de vida).
- Especiação em duas etapas: surgimento contínuo de espécies incipientes e persistência delas por
  um período de transição; estrutura espacial e dispersão regulam as duas.

## Januario, Pinsky & Rabosky (2025)

- Arcabouço para testar se a persistência metapopulacional (ocupação de equilíbrio e
  resiliência) prediz a taxa de especiação nas pontas.
- Peixes demersais: um índice de persistência tem sinal filogenético, mas nenhum prediz a taxa
  de especiação; os autores falam em desacoplamento entre a dinâmica de décadas e a de milhões
  de anos.
- Os dois trabalhos tratam de taxas de diversificação; aqui as taxas são de evolução de traits.

## Leibold, Govaert, Loeuille, De Meester & Urban (2022)

- A evolução adaptativa altera a montagem local e regional, e o efeito depende da
  heterogeneidade ambiental e da sua estrutura espacial. Resultado típico: efeitos de prioridade
  mediados pela evolução (colonizadores precoces se adaptam e alteram o estabelecimento dos que
  chegam depois; monopolização da comunidade).
- Três taxas relativas decidem o desfecho (Fig. 1): colonização, regulação da comunidade e
  evolução. Evolução lenta: quadro ecológico clássico. Evolução rápida em relação à colonização:
  monopolização e, no extremo, radiação adaptativa. Colonização alta: fluxo gênico afoga a
  adaptação local (gene swamping).
- Competição: deslocamento de caráter, convergência que produz coocorrência quase neutra quando
  há mais espécies que nichos, ciclos eco-evolutivos. Interações tróficas: perda de defesas na
  ausência do predador (genetic subsidy).
- Lacuna apontada pelos autores: os modelos macroecológicos não incorporam adequadamente quanto
  a adaptação local e a montagem alteram as taxas de especiação e extinção e, por elas, as
  filogenias e a biota regional.
- Cenários sugeridos: frequência de extinção de manchas (distúrbios reiniciam a montagem; em
  excesso impedem a adaptação), autocorrelação espacial do ambiente, ordem de chegada,
  adaptações secundárias que revertem a monopolização.
- Métodos de partição de variância de traits: estatísticas T de Violle et al. (2012), Lepš et
  al. (2011). As estatísticas T estão em `cati::Tstats()` (CRAN 0.99.6, GPL >= 2; pode ir para
  `Suggests`): T_IP.IC (variância dentro das populações sobre a variância dentro das
  comunidades), T_IC.IR (comunidades sobre o pool regional) e T_PC.PR (entre comunidades sobre o
  pool regional), cada uma com seu modelo nulo (local, regional.ind, regional.pop). Exigem traits
  de indivíduos, com fatores de parcela e de espécie.
- Interpretação nossa: eles falam de taxas de evolução relativas à colonização em escala
  ecológica, não de taxas por ramo numa filogenia.

## Leibold, Urban, De Meester, Klausmeier & Vanoverbeke (2019)

- Pergunta: o papel da evolução adaptativa local na coexistência, numa paisagem de manchas
  conectadas sujeitas a distúrbio.
- Dois modelos. Um modelo de ocupação de manchas (determinístico, sem coexistência local; o
  residente adaptado não pode ser invadido), resolvido no Mathematica. E um modelo baseado em
  indivíduos numa grade 10 × 10 com dois tipos de mancha distribuídos ao acaso (ótimos 0,4 e
  0,6), trait de 10 locos dialélicos com mutação, aptidão gaussiana, extinção de manchas com
  probabilidade Pext e dispersão com núcleo exponencial (código em Java no Zenodo,
  10.5281/zenodo.2280982, conforme o artigo).
- Resultado: os colonizadores precoces se adaptam à mancha antes da chegada dos pré-adaptados
  (efeitos de prioridade mediados pela evolução). A repetição disso faz as duas espécies
  convergirem para distribuições bimodais de trait iguais na escala regional e para ocupações
  iguais dos dois tipos de mancha, embora cada mancha continue monopolizada por uma espécie.
  Os autores chamam isso de neutralidade regional: nicho na escala local, equivalência na
  regional, e deriva regional nas abundâncias relativas.
- Mais neutralidade regional com dispersão baixa (menos de 100 indivíduos dispersando por
  mancha por geração, Fig. 5); com extinção baixa e dispersão alta ela demora mais.
- Partição de variância: a neutralidade regional reduz a fração ambiental da composição e
  infla a espacial e o resíduo. Para os traits (e não a composição), o ambiente explica quase
  toda a variação. A partição da composição, sozinha, não separa neutralidade puramente
  ecológica de neutralidade regional evoluída.
- Previsão: os táxons mostram pouco rastreamento ambiental, mas os valores de trait locais
  (independentes do táxon) mostram muito. Isso exige medir traits no nível da população.

## Urban et al. (2020)

- Síntese de 500 estudos eco-evolutivos com 14 mecanismos pelos quais a evolução local
  amortece, amplifica ou cria padrões ecológicos espaciais.
- A evolução amortece os padrões em 85% dos estudos, amplifica em 13% e cria em 8% (a soma
  passa de 100% porque um estudo pode ter mais de um padrão). A
  adaptação local iguala a aptidão entre ambientes e alisa a heterogeneidade.
- Montagem: colonizadores precoces podem se adaptar e alterar permanentemente riqueza e
  composição. A evolução de neutralidade é citada como mecanismo pouco visto empiricamente.
- Questão em aberto 5: em que condições a especiação converte diferenças entre populações em
  diferenças ecológicas entre espécies. A questão 2 lembra que o resultado depende das taxas
  relativas de evolução, colonização e redundância de nicho.

## Ponte Leibold–Januario (interpretação nossa)

- Leibold et al. (2022) dizem que falta incorporar à macroecologia como a montagem e a
  adaptação local alteram especiação e extinção. Januario et al. (2025) testaram uma ligação
  desse tipo (persistência metapopulacional → especiação nas pontas) e acharam desacoplamento.
- Leibold et al. (2019) e Urban et al. (2020) oferecem mecanismos para esse desacoplamento:
  - a neutralidade regional deixa as espécies equivalentes na escala regional, de modo que
    propriedades da dinâmica local, como a persistência, não se traduzem em diferenças entre
    espécies;
  - a evolução local amortece os padrões espaciais;
  - a divergência adaptativa fica entre populações, e só vira diferença entre espécies se
    houver especiação (Urban et al., questão 5).
- Para taxas de evolução de traits, a previsão análoga é: com evolução local rápida em
  relação à colonização e dispersão baixa, a evolução entre populações é rápida, mas as médias
  das espécies convergem. As taxas por ramo estimadas com médias de espécies ficam baixas e
  pouco ligadas à dinâmica local; as estatísticas T (T_IP.IC, T_PC.PR) carregam o sinal.
- O que o comdynamite pode testar com verdade conhecida: variar extinção de manchas e
  dispersão (os eixos da Fig. 5 de Leibold et al. 2019) e medir se a história de colonização
  (`history_rates()`) deixa assinatura nas taxas estimadas por ramo, comparando traits de
  espécies e de populações. Nenhum dos três trabalhos usa taxas por ramo; a ligação é nossa.
- Limite: o comdynamite ainda não simula especiação acoplada à metacomunidade, então a
  conversão de divergência entre populações em especiação (o elo de Januario) fica fora por
  enquanto. Plano para fechar esse limite em `dev/plano-especiacao.md`.

### Experimento da ponte (resultados, `data-raw/bridge_experiments.R`)

Árvore de 30 espécies, grade 10 × 10, 600 passos; dispersão {0,001; 0,01; 0,1} × extinção de
manchas {0; 0,01; 0,05} × autocorrelação do ambiente (`env_bandwidth` {0,5; 3}), 20 réplicas
(360 simulações); h {0; 0,5} × c {0; 1; 3}; desenhos de 20 ou 100 sítios × 5 ou 20 indivíduos
por população; `intraspecific_sd = 5`, `population_sd = 1`. Bootstrap (nboot = 99) num
subconjunto. 48 min em 6 núcleos.

- Dispersão 0,001 com extinção 0,05 colapsa a metacomunidade (sem indivíduos amostrados): 40
  simulações sem dados na camada de indivíduos. É o limite em que a extinção vence a
  recolonização.
- Taxas por ramo: a correlação de Spearman entre a taxa RRphylo das pontas e o turnover de
  parceiros fica entre −0,01 e 0,20 (média por célula), maior com dispersão intermediária ou
  alta. O rastro da montagem nas taxas por ramo, mesmo plantado por `turnover_rates()`, é
  fraco, como na vignette. (A correlação em todos os ramos não é interpretável:
  `turnover_rates(window = 0.1)` só muda as pontas e os ramos internos ficam todos com taxa 1.)
- h é recuperado em todos os desenhos (RMSE 0,01–0,04).
- c é identificável nestes desenhos, ao contrário de SJC: poder para c = 1 (contra o nulo c = 0
  da mesma célula) de 0,78 com 20 sítios × 5 indivíduos e w livre, e ≥ 0,94 nos outros
  desenhos; com w fixo, ≥ 0,997. Os desenhos com 20 sítios × 5 indivíduos perdem poder em
  metacomunidades com dispersão muito baixa ou extinção (0,53–0,68 em várias células). O erro
  de ĉ cai sobretudo com mais indivíduos por população (5 → 20) e com mais espécies por sítio;
  a beta diversidade da composição não importa.
- Com w livre, c forte é atenuado (c = 3: mediana 1,85–2,43); com w fixo quase não (2,81–2,97).
- O aviso de w mal identificado aparece em 95–99% das réplicas sem deslocamento e quase nunca
  com c ≥ 1: funciona como diagnóstico de ausência de deslocamento.
- Estatísticas T (nulo por espécie, 100 sítios × 20 indivíduos): SES de T_IP.IC 0,08 sem
  processo, 4,15 só com adaptação local, −0,8 a −1,0 só com deslocamento e 0,25 com h = 0,5 e
  c = 3: a compensação entre os dois processos aparece também com o nulo por espécie.
- Bootstrap de `population_shifts()` (versão do PR #57): tamanho 0,025 com 20 sítios × 5
  indivíduos, mas 0,29 com 100 sítios × 20 indivíduos; cobertura dos intervalos de c = 3 de
  0,04–0,13.
- **Causa e correção (PR #58).** O gerador nulo usava os slopes do ajuste com deslocamento; o
  erro deles (ĥ 0,548 contra 0,5) ia para os resíduos e inflava o dp de população (1,57 contra
  1). Desvios de população viesam ĉ para cima quando o kernel é largo (a derivada do push em
  relação à própria média é exp(−g²/2w²)(1 − g²/w²), positiva para g < w), então o nulo
  exagerava o viés (média das réplicas nulas 0,21–0,30 contra 0,13–0,14 na distribuição
  verdadeira de ĉ numa mesma metacomunidade). Correção: nulo a partir do ajuste sem
  deslocamento, constantes das populações que entram no ajuste, variância de população por
  momentos ponderados pela precisão e, no conjunto ajustado, variância de amostragem escalada por
  (1 − c D)², já que o erro da média entra no push removido (achado da auditoria). A
  documentação dizia que o viés por desvios de população era para baixo; o sinal depende de w.
- Nova rodada do subconjunto (`data-raw/bridge_bootstrap.R`, novas simulações, 81 por
  desenho, 28 min):

  | c | desenho | rejeição a 5% (antes) | cobertura de c (antes) | cobertura de h (antes) | ĉ corrigido mediano (antes) |
  |---|---|---|---|---|---|
  | 0 | 20 × 5 | 0,06 (0,03) | 0,90 (0,91) | 0,87 (0,89) | −0,02 (−0,04) |
  | 0 | 100 × 20 | 0,00 (0,29) | 0,68 (0,80) | 0,64 (0,87) | 0,06 (−0,01) |
  | 3 | 20 × 5 | 0,93 (0,99) | 0,41 (0,13) | 0,88 (0,04) | 2,34 (3,02) |
  | 3 | 100 × 20 | 0,99 (1,00) | 0,49 (0,04) | 0,78 (0,02) | 2,72 (3,36) |

  O teste deixou de rejeitar demais; em 100 × 20 pareceu conservador (0 de 81 rejeições,
  mediana do p 0,58). Diagnóstico à parte (100 × 20, extinção 0,01, 6 metacomunidades): a forma
  do nulo está certa (150 redraws verdadeiros contra o nulo de um conjunto de dados: rejeição
  0,03–0,11, média 0,055), mas cada conjunto de dados contra o próprio nulo rejeita 3,3% (180
  conjuntos; somando a rodada acima, 6 de 261 ≈ 2,3%). O centro do nulo acompanha a estimativa
  (correlação entre ĉ e a média das réplicas nulas 0,2–0,33), sobretudo pelo dp de população,
  estimado dos mesmos resíduos que movem ĉ. Decisão: documentar o leve conservadorismo (PR
  #58) em vez de calibrar por bootstrap duplo (custo nboot²). Os intervalos de c = 3 cobrem mais, mas ainda bem abaixo de 95%, e a correção
  agora fica curta em vez de passar do ponto: o viés de ĉ muda com c, como a documentação
  avisa. Sob c = 0 em 100 × 20, ĥ tem viés +0,03 com intervalo de largura 0,05, e a cobertura
  de c e h cai porque o conjunto ajustado é gerado com c0 = ĉ > 0. Recomendação: usar
  `p_displacement` como guia; tratar intervalos e correção de c como aproximados.

#### Por que c é identificável na ponte e não em SJC (`data-raw/bridge_narrow.R`)

- Mesmo número de sítios (20), mas o desvio-padrão de ĉ sob o nulo é 1,38 na calibração de SJC
  e 0,32 no desenho 20 × 5 da ponte.
- O que governa a identificabilidade é a largura do kernel em relação ao espaçamento das espécies
  no trait dentro dos sítios (`spread`, mediana do desvio-padrão das médias de população por
  sítio), não o ruído de medida nem o tamanho da amostra. Descritores adimensionais (uma réplica
  da ponte com h = 0,5 e c = 0; em SJC, w são os valores da calibração, não estimativas):

  | desenho | pops usadas | espécies/sítio | ruído de meia média / spread | w / spread | dp do push líquido do ambiente / spread (por unidade de c) | EP aproximado de ĉ (w fixo) |
  |---|---|---|---|---|---|---|
  | SJC, w = 0,1 | 88 | 6,1 | 0,11 | 0,20 | 0,014 | 0,58 |
  | SJC, w = 0,25 | 88 | 6,1 | 0,11 | 0,50 | 0,056 | 0,15 |
  | ponte 20 × 5 | 181 | 10,7 | 0,35 | 1,83 | 0,227 | 0,08 |
  | ponte 100 × 20 | 923 | 10,9 | 0,16 | 1,81 | 0,192 | 0,02 |

  O ruído relativo é até menor em SJC. A diferença está no push: com kernel estreito em relação ao
  espaçamento, poucos pares de populações ficam a menos de w, o push quase não varia entre as
  populações de uma espécie e c e w se compensam.
- Teste causal: o mesmo desenho 20 × 5 da ponte, com o kernel encolhido. Poder (contra o
  quantil 95% do nulo da mesma razão; 60 réplicas):

  | w / spread | c = 1, w livre | c = 1, w fixo | c = 3, w livre | c = 3, w fixo |
  |---|---|---|---|---|
  | 0,2 | 0,10 | 0,08 | 0,18 | 0,17 |
  | 0,5 | 0,33 | 0,53 | 0,75 | 0,98 |
  | 1,8 | 0,98 | 1,00 | 1,00 | 1,00 |

- A ponte planta w = média de `sigma` (a largura de nicho), larga em relação ao espaçamento dos
  ótimos dentro dos sítios: um deslocamento forte e fácil de ver. Isso é escolha do simulador,
  não propriedade do método.
- Em SJC o melhor w fica no piso da grade (0,10 ≈ 0,2 × spread) e o perfil é plano. Os dados são
  compatíveis com ausência de deslocamento ou com deslocamento por um kernel estreito, e o desenho
  não separa as duas coisas. Na ponte, sob c = 0, o aviso de w mal identificado também dispara em
  95–99% das réplicas: w não identificado é o sintoma esperado de c ≈ 0, não uma falha específica
  de SJC.
- Interpretação nossa: em anuros adultos o CRC difere entre espécies co-ocorrentes muito
  mais do que a largura plausível de um kernel de competição por tamanho, então w / spread
  pequeno é o caso provável em dados reais desse tipo. Para identificar c, ajuda (i) fixar w por
  conhecimento independente, (ii) usar traits em que as espécies co-ocorrentes se sobrepõem, ou
  (iii) desenhos com mais espécies próximas no trait por sítio.

- Ressalva de escala: h e c entram na amostragem, não evoluem na dinâmica; as taxas por ramo
  vêm de `turnover_rates()`; as duas camadas não compartilham o trait. O experimento testa a
  recuperação de sinais plantados sob estruturas de montagem realistas, não se a montagem causa
  evolução, e o desacoplamento de Januario et al. (2025) fica fora (proposta de teste de poder
  em separado).

## McPeek (2017), *Evolutionary Community Ecology*

Leitura feita por subagente; três citações conferidas no texto extraído.

- Mecanismo: a mudança do traço médio é V_z vezes a soma, sobre as espécies que interagem, da
  abundância de cada uma vezes o gradiente de seleção (eq. 3.7). V_z é constante e pequeno; toda
  variação de velocidade vem do contexto da comunidade (quem está presente, abundâncias, traits).
- Não usa o vocabulário de métodos comparativos (BM, sigma², early burst, taxas por ramo). A
  ligação com modelos de taxas por ramo é interpretação nossa: o livro dá o mecanismo, não a
  estatística macroevolutiva.
- Evolução acelera: invasão de uma comunidade (mudança rápida que desacelera até o equilíbrio,
  Fig. 3.2–3.3), mudanças de hábitat ou hospedeiro ("must occur very rapidly, probably over, at
  most, tens of generations", p. 149), peak shifts na diferenciação, ciclos de traits,
  reembaralhamento das comunidades no Quaternário. "Large and changing types of selection should
  be expected most when species are first adapting to a new selection regime" (p. 129).
- Evolução desacelera: equilíbrios evolutivos (estase com seleção direcional nos componentes),
  espécies neutras (até 12 *Enallagma* por cerca de 10 milhões de anos), adição de predador
  (enfraquece a seleção sobre o recurso), fluxo gênico.
- Deslocamento de caráter não é inevitável: exige abundâncias altas, alta produtividade,
  interações fortes; propaga-se de recursos para consumidores e predadores.
- Metacomunidade (cap. 6): dispersão homogeneiza abundâncias e traits; acima de um limiar de
  diferença entre ótimos locais há transição abrupta para baixa dispersão e diferenciação local.
  Os modelos do cap. 6 têm uma espécie em duas ou J manchas, não uma metacomunidade
  multiespecífica com traits evoluindo.
- "By reconstructing trait evolution across multiple interacting clades, hypotheses about
  coevolutionary responses among species can also be tested" (p. 274).

## Cenários de simulação derivados (interpretação nossa)

1. Explosão na colonização: mudança de taxa no ramo em que a linhagem coloniza um novo tipo de
   mancha, decaindo ao longo do ramo; a linhagem residente responde com atraso.
2. Mudança de hábitat como gatilho: classe de taxa rápida em `trait_mk()` (ou mudança em
   `trait_rates()`) nos ramos com transição de hábitat, seguida de estase.
3. Taxa como função do ambiente: taxa por linhagem crescendo com a produtividade ou a força das
   interações nas manchas ocupadas.
4. Posição trófica: recursos mais rápidos que consumidores, consumidores mais rápidos que
   predadores; queda da taxa do recurso quando um predador entra.
5. Estase pontuada por mudanças na composição: taxa do ramo proporcional ao número de entradas e
   saídas de parceiros de interação ao longo dele, contra um BM homogêneo como nulo. É a previsão
   mais defensável de que a montagem deixa rastro nas taxas.
6. Taxas acopladas entre clados que interagem (recurso e consumidor).
7. Dispersão: linhagens generalistas de alta dispersão com taxas baixas; especialistas de baixa
   dispersão divergem; variar a diferença entre ótimos para achar o limiar.
8. Mudanças de taxa por Poisson com intensidade variável no tempo, concentradas em eventos de
   reembaralhamento das comunidades e compartilhadas entre linhagens simpátricas, contra mudanças
   independentes entre linhagens.
9. Frequência de extinção de manchas e autocorrelação espacial do ambiente (Leibold et al. 2022),
   cruzadas com a dispersão como na Fig. 5 de Leibold et al. (2019), para achar onde a
   neutralidade regional apaga o rastro da montagem nas taxas entre espécies.
10. Variação intraespecífica, adaptação local e deslocamento de caráter. Indivíduos gerados a
    partir das abundâncias simuladas com `sample_individuals()`, com traits ajustados ao ambiente
    do sítio (`local_adaptation`, h) e afastados dos co-ocorrentes (`displacement`, c), e
    estatísticas T com `tstats()` (`data-raw/demo_tstats.R`).
    - Previsão corrigida. A adaptação local *sobe* T_IP.IC: as populações convergem para o
      ambiente do sítio, as espécies co-ocorrentes ficam parecidas e a variância dentro da
      comunidade encolhe (0,02 → 1,0), enquanto T_IC.IR e T_PC.PR caem (0,54 e 0,89 → < 0,01).
      Só o deslocamento de caráter *baixa* T_IP.IC (e sobe T_PC.PR). A previsão original
      ("adaptação local forte deve baixar T_IP.IC") estava errada.
    - Os dois processos se compensam: com h = 0,9, `displacement = 10` leva T_IP.IC de 0,61 a
      0,011 e T_PC.PR de 0,010 a 0,59, e "nem adaptação nem deslocamento" fica parecido com "os
      dois fortes" (T_PC.PR 0,89 contra 0,59; T_IP.IC 0,017 contra 0,011). As razões sozinhas não
      separam os processos.
    - T_IP.IC só deve ser reportado com `null_ip = "species"`, que troca as médias das
      populações entre populações da mesma espécie e mantém composição, abundâncias e variâncias
      dentro das populações. O nulo por sítio (permutar indivíduos dentro do sítio, Taudiere &
      Violle 2016) só pergunta se as populações diferem: dá efeito negativo grande sempre que
      co-ocorrentes têm médias diferentes, mesmo sem nenhuma resposta ao sítio. Com o nulo por
      espécie, o efeito padronizado é positivo sob convergência (adaptação local a um ambiente
      comum) e negativo sob divergência (deslocamento). Buscas no Consensus, Scholar Gateway e OpenAlex não acharam
      esse nulo na literatura das estatísticas T; não é prova de que não exista
      [CITAR: nulo restrito por espécie para T_IP.IC, se houver].
    - A inferência sobre os processos vem de `population_shifts()`, que ajusta h e c juntos às
      médias das populações, dentro de espécie, com o ambiente do sítio e a abundância dos
      competidores. No piloto ABC (`data-raw/pilot_tstats_identifiability.R`), as três razões
      recuperam h (R² 0,87–0,89) mas mal c (0,56–0,58), e a regressão conjunta recupera os dois
      (0,99 e 0,95–0,96); o piloto usa o mesmo modelo que gera os dados, então é um limite
      superior. A abordagem segue Goldberg & Lande (2006, Evolution 60:1344–1357,
      doi:10.1554/05-696.1) e Adams & Collyer (2007, Evolution 61:510–515,
      doi:10.1111/j.1558-5646.2007.00063.x), que pedem avaliar a divergência entre simpátricas
      depois de descontar o gradiente ambiental (verificadas no Crossref).
    - Limite em dados reais: no desenho de São José dos Campos (20 poças, 15 espécies, ≤ 5
      indivíduos por população; calibração por simulação do próprio desenho), h é recuperado sem
      viés (±0,01 por desvio-padrão do gradiente), mas c não é identificável com a largura do
      kernel livre: sem deslocamento, ĉ varia de −1,9 a 2,7 (90% central), e o poder para c = 4
      fica em 30–36%. Com a largura fixa no valor verdadeiro, c é bem estimado; o problema é
      identificar w. Desvios de populações inteiras também enviesam ĉ para baixo, o que o
      bootstrap de `population_shifts()` incorpora no p-valor centrado.
    - Para o artigo: as estatísticas T precisam do ambiente e da composição local para separar
      os processos, e o experimento da ponte (abaixo) deve mapear em que desenhos amostrais e em
      que metacomunidades (dispersão, extinção de manchas, autocorrelação do ambiente) c passa a
      ser estimável.
    Interpretação nossa; não verifiquei se a literatura das estatísticas T já trata da
    compensação entre adaptação local e deslocamento
    [CITAR: trabalhos sobre T_IP/IC e deslocamento de caráter ou plasticidade].

Para os cenários 1, 2, 5 e 8 o pacote já tem quase tudo: `simulate_metacommunity()` dá a
história de ocupação e `trait_rates(shift_edges = )` ou `trait_mk(edge_rates = )` aceitam taxas
definidas por ramo, e `history_rates()` e `turnover_rates()` traduzem a história da metacomunidade
em taxas por ramo. Os resultados estão no artigo `vignettes/articles/rates-assembly.Rmd`.
