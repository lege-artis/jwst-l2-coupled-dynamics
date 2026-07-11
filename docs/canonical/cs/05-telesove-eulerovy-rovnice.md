# Kapitola 5 — Tělesové Eulerovy rovnice jako specializace

> **Předpokládané čtení.** Kapitola 02 (kinetická energie na $SE(3) \times SE(3)$, zejména §2.7 levá invariance rotační metriky a §2.7.1 explicitní konstrukce). Kapitola 03 (moment síly z gravitačního gradientu, zejména §3.7 zarámovaný tělesový výraz a §3.7.2 orientační Jacobiho matice + tabulka stability librační rovnováhy po opravě znaménka OQ-FORT-1). Kapitola 04 (Euler-Lagrange + Euler-Poincarého redukce; tato kapitola je *specializací* zarámovaného výsledku kapitoly 04 §4.4 na soustavu hlavních os, s explicitním momentem síly z gravitačního gradientu dosazeným z kap. 03 §3.7).
>
> **Stav.** v0.1 (kanonická úroveň, 3. kolo, 2026-05-25). Všechny čtyři vnitřní otevřené otázky kapitoly OQ-5.1, OQ-5.2, OQ-5.3, OQ-5.4 vyřešeny v tomto kole; viz §5.8 pro souhrn uzávěru.

> **Poznámka ke způsobu a předpokladům odvození.** Tato CS kanonická kapitola (v0.1) je odvozena z anglické kanonické kapitoly `docs/canonical/en/05-body-frame-euler-equations.md` (v0.1, 3. kolo). EN verze zůstává autoritativním pramenem při budoucích aktualizacích.

---

## §5.1 Co tato kapitola dělá

Kapitola 04 odvodila, skrze Euler-Poincarého redukci levo-invariantního rotačního Lagrangiánu na každém faktoru $SO(3)$, obecné tělesové pohybové rovnice

$$
\mathbf{I}_X \dot{\boldsymbol{\Omega}}_X \;+\; \boldsymbol{\Omega}_X \times (\mathbf{I}_X \boldsymbol{\Omega}_X) \;=\; \boldsymbol{\tau}_X^{\mathrm{body}}, \qquad X \in \{A, B\},
$$

kde $\boldsymbol{\Omega}_X$ je tělesová úhlová rychlost tělesa $X$, $\mathbf{I}_X$ je jeho (v tělesové soustavě konstantní) tenzor setrvačnosti a $\boldsymbol{\tau}_X^{\mathrm{body}}$ je tělesový vnější moment síly na těleso $X$. **Tato kapitola specializuje onen obecný výsledek na soustavu hlavních os** (kde je $\mathbf{I}_X$ diagonální) a dosazuje explicitní moment síly z gravitačního gradientu z kapitoly 03 §3.7 k vyprodukování uzavřených skalárních pohybových rovnic, jež implementace inženýrské úrovně propaguje a fortranový lege-artis backend validuje.

Kapitola rovněž vynáší **specifické strukturální rysy tvaru hlavních os**, jež jsou důležité pro numerickou implementaci: antisymetrickou vazební strukturu Eulerových členů $(I_2 - I_3)\Omega_y\Omega_z$, zvláštní roli sférické a axisymetrické třídy symetrie (kde některé Eulerovy členy identicky mizí) a spojení s testovacími případy analytických limitů (bezsilový symetrický setrvačník, Lagrangeův setrvačník, nestabilita asymetrického setrvačníku Džanibekova), vůči nimž bude kapitola 06 validovat.

Tato kapitola sleduje lege-artis metodologii stanovenou v kapitolách 02–04: každá rovnice je odvozena z výsledků předchozí kapitoly explicitně, každá specializace třídy symetrie je zdůvodněna, každý analytický limit je pojmenován a ukázán na svůj validační testovací případ kapitoly 06. Kapitola je strukturována jako **uzavřená specializační reference** — čtenář, jenž zná výsledek kapitoly 04 (tělesovou Eulerovu rovnici v souřadnicově nezávislé formě), může přijít k této kapitole pro samotný skalární tvar hlavních os.

## §5.2 Specializace tělesových Eulerových rovnic na hlavní osy

Zarámovaná tělesová Eulerova rovnice kap. 04 §4.4 zní, v souřadnicově nezávislé formě:

$$
\mathbf{I}_X \dot{\boldsymbol{\Omega}}_X \;+\; \boldsymbol{\Omega}_X \times (\mathbf{I}_X \boldsymbol{\Omega}_X) \;=\; \boldsymbol{\tau}_X^{\mathrm{body}}.
$$

Tenzor setrvačnosti $\mathbf{I}_X$ je symetrický a pozitivně definitní (kap. 02 §2.4), takže je diagonalizovatelný ortogonální změnou báze. **Soustava hlavních os** je tělesová ortonormální soustava $(\hat{\mathbf{e}}_1^X, \hat{\mathbf{e}}_2^X, \hat{\mathbf{e}}_3^X)$, v níž $\mathbf{I}_X$ nabývá diagonálního tvaru $\mathbf{I}_X = \mathrm{diag}(I_1^X, I_2^X, I_3^X)$, s **hlavními momenty setrvačnosti** $I_1^X, I_2^X, I_3^X$ jako diagonálními prvky (vlastními hodnotami $\mathbf{I}_X$). Rotace do soustavy hlavních os je podobnostní transformací aplikovanou na tělesovou rovnici; v této soustavě pracujeme napříč touto kapitolou.

### §5.2.1 Rozvoj složka po složce (řeší OQ-5.1)

Zapiš zarámovanou rovnici po složkách v soustavě hlavních os. Nechť $\boldsymbol{\Omega}_X = (\Omega_x, \Omega_y, \Omega_z)$ v bázi hlavních os (potlačujíc index $X$ na jednotlivých složkách pro čitelnost; $X$ je obnoveno při nejednoznačnosti). Dva členy vlevo jsou:

**Člen 1 — tělesové úhlové zrychlení:**

$$
\mathbf{I}_X \dot{\boldsymbol{\Omega}}_X
\;=\;
\mathrm{diag}(I_1, I_2, I_3)
\begin{pmatrix} \dot\Omega_x \\ \dot\Omega_y \\ \dot\Omega_z \end{pmatrix}
\;=\;
\begin{pmatrix} I_1 \dot\Omega_x \\ I_2 \dot\Omega_y \\ I_3 \dot\Omega_z \end{pmatrix}.
$$

**Člen 2 — vektorový součin úhlová rychlost / moment hybnosti:** $\mathbf{I}_X \boldsymbol{\Omega}_X = (I_1 \Omega_x, I_2 \Omega_y, I_3 \Omega_z)$, takže:

$$
\boldsymbol{\Omega}_X \times (\mathbf{I}_X \boldsymbol{\Omega}_X)
\;=\;
\begin{vmatrix} \hat{\mathbf{e}}_1 & \hat{\mathbf{e}}_2 & \hat{\mathbf{e}}_3 \\ \Omega_x & \Omega_y & \Omega_z \\ I_1\Omega_x & I_2\Omega_y & I_3\Omega_z \end{vmatrix}
\;=\;
\begin{pmatrix}
\Omega_y(I_3\Omega_z) - \Omega_z(I_2\Omega_y) \\
\Omega_z(I_1\Omega_x) - \Omega_x(I_3\Omega_z) \\
\Omega_x(I_2\Omega_y) - \Omega_y(I_1\Omega_x)
\end{pmatrix}
\;=\;
\begin{pmatrix}
(I_3 - I_2)\,\Omega_y\Omega_z \\
(I_1 - I_3)\,\Omega_z\Omega_x \\
(I_2 - I_1)\,\Omega_x\Omega_y
\end{pmatrix}.
$$

Sečtením dvou členů a přirovnáním k tělesovému momentu síly $\boldsymbol{\tau}_X^{\mathrm{body}} = (\tau_x^{\mathrm{body}}, \tau_y^{\mathrm{body}}, \tau_z^{\mathrm{body}})$:

$$
\boxed{
\begin{aligned}
I_1 \dot\Omega_x \;-\; (I_2 - I_3)\, \Omega_y \Omega_z \;&=\; \tau_x^{\mathrm{body}}, \\
I_2 \dot\Omega_y \;-\; (I_3 - I_1)\, \Omega_z \Omega_x \;&=\; \tau_y^{\mathrm{body}}, \\
I_3 \dot\Omega_z \;-\; (I_1 - I_2)\, \Omega_x \Omega_y \;&=\; \tau_z^{\mathrm{body}}.
\end{aligned}
}
$$

Toto jsou klasické **Eulerovy dynamické rovnice z roku 1758** pro tuhé těleso v jeho soustavě hlavních os, obnovené zde jako strukturální specializace sady rovnic z kapitoly 04 odvozené Euler-Poincarého redukcí spíše než jako Newtonovsko-Eulerovo tvrzení o rovnováze sil. Obě cesty souhlasí (křížová validace kapitoly 04 §4.5); tato kapitola bere sadu rovnic odvozenou Euler-Poincaré jako kanonické tvrzení a dokazuje, že se specializuje na týž algebraický tvar, jejž Euler zapsal roku 1758. **(Řeší OQ-5.1.)**

Rovnice tělesa-$B$ jsou symetrické s $A \leftrightarrow B$ a nezávislou sadou tělesových hlavních momentů a složek úhlové rychlosti. Dvě sady tělesových rovnic jsou **vázány pouze skrze orientačně závislé momenty sil z gravitačního gradientu** $\boldsymbol{\tau}_X^{\mathrm{body}}$, jež závisí na relativní orientaci obou těles skrze tělesový směr $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, X}$ (kapitola 03 §3.7 zarámovaný tvar).

### §5.2.2 Kontrola zdravého rozumu konvence znaménka

Člen vektorového součinu vstupuje s **minusovým** znaménkem na levé straně zarámovaných rovnic (ekvivalentně s **plusovým** znaménkem na pravé straně, je-li přesunut). Toto znaménko pochází z odvození kap. 04 §4.4 $\mathrm{ad}^*_{\boldsymbol{\Omega}}\,\boldsymbol{\mu} = -\boldsymbol{\Omega} \times \boldsymbol{\mu}$, jež se propaguje skrze Euler-Poincarého rovnici $d(\mathbf{I}\boldsymbol{\Omega})/dt = \mathrm{ad}^*_{\boldsymbol{\Omega}}(\mathbf{I}\boldsymbol{\Omega}) + \boldsymbol{\tau}^{\mathrm{body}}$ a dává $\mathbf{I}\dot{\boldsymbol{\Omega}} = -\boldsymbol{\Omega} \times (\mathbf{I}\boldsymbol{\Omega}) + \boldsymbol{\tau}^{\mathrm{body}}$ nebo ekvivalentně $\mathbf{I}\dot{\boldsymbol{\Omega}} + \boldsymbol{\Omega} \times (\mathbf{I}\boldsymbol{\Omega}) = \boldsymbol{\tau}^{\mathrm{body}}$.

**Křížová kontrola znaménkem nestability Džanibekova.** Pro asymetrické těleso s $I_1 < I_2 < I_3$ rotující převážně kolem prostřední osy $\hat{\mathbf{e}}_2$ — tj. $\boldsymbol{\Omega} = (0, \Omega_2^0, 0) + \boldsymbol{\delta\Omega}$ s $|\boldsymbol{\delta\Omega}| \ll |\Omega_2^0|$ — dávají linearizované rovnice pro $\delta\Omega_x, \delta\Omega_z$ odvozené dosazením do zarámovaných §5.2.1 rovnic a ponecháním členů prvního řádu vlastní hodnotu nestability $\sigma^2 = (\Omega_2^0)^2 (I_2 - I_1)(I_3 - I_2)/(I_1 I_3) > 0$ (kladnou, tudíž nestabilní). Pro rotaci kolem kterékoli krajní osy ($\hat{\mathbf{e}}_1$ s $\boldsymbol{\Omega} = (\Omega_1^0, 0, 0)$ + perturbace nebo $\hat{\mathbf{e}}_3$ podobně) je analogická vlastní hodnota $\sigma^2 = -(\Omega_1^0)^2 (I_2 - I_1)(I_3 - I_1)/(I_2 I_3) < 0$ (záporná, tudíž stabilní). Znaménko zarámovaných Eulerových rovnic je tedy konzistentní s učebnicovou klasifikací nestability tenisové rakety (Landau-Lifšic, sv. I §37; Goldstein 2002 §5.6); testovací případ Džanibekova (kap. 06 §6.5 / `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md` §2.2) je trvalou regresní pojistkou.

## §5.3 Specializace tříd symetrie

Struktura Eulerových rovnic se dramaticky zjednodušuje, má-li těleso dodatečnou symetrii. Tři relevantní případy:

| Třída symetrie tělesa | Hlavní momenty | Struktura Eulerovy rovnice |
|---|---|---|
| **Sférická** | $I_1 = I_2 = I_3 \equiv I_0$ | Všechny Eulerovy smíšené členy identicky mizí; rovnice se redukují na $I_0 \dot{\boldsymbol{\Omega}} = \boldsymbol{\tau}^{\mathrm{body}}$ — čistě momentem-síly-hnané úhlové zrychlení bez setrvačné vazby. Bezsilový případ ($\boldsymbol{\tau} = 0$) dává $\boldsymbol{\Omega} = \mathrm{const}$. Sférická tělesa nemají *žádnou vnitřní rotační dynamiku* v tělesové soustavě — veškerá složitost je ve vývoji orientace. |
| **Axisymetrická** ($I_1 = I_2 \equiv I_\perp$, $I_3 \equiv I_\parallel$) | Dva odlišné momenty | Rovnice $\dot\Omega_z$ se rozpojí: $I_\parallel \dot\Omega_z = \tau_z^{\mathrm{body}}$. Dvojice $(x, y)$ se váže jako $I_\perp \dot\Omega_x = (I_\perp - I_\parallel) \Omega_y \Omega_z + \tau_x^{\mathrm{body}}$ (a cyklicky). Bezsilový případ je **volný symetrický setrvačník** s regulární precesí $\boldsymbol{\Omega}$ kolem osy symetrie rychlostí $\lambda = (I_\parallel - I_\perp)\Omega_z / I_\perp$. |
| **Asymetrická** (všechny tři $I_i$ odlišné) | Tři odlišné momenty | Všechny tři Eulerovy smíšené členy jsou nenulové. Bezsilový případ vykazuje **efekt tenisové rakety Džanibekova**: rotace kolem prostřední hlavní osy ($I_2$ mezi $I_1$ a $I_3$) je nestabilní; malé perturbace narůstají v periodické 180° přeskoky osy spinu. |

### §5.3.1 Sférický případ — vypracované detaily

Ve sférickém případě $I_1 = I_2 = I_3 = I_0$ každý koeficient $(I_i - I_j)$ mizí, člen vektorového součinu $\boldsymbol{\Omega} \times (\mathbf{I}\boldsymbol{\Omega}) = \boldsymbol{\Omega} \times (I_0 \boldsymbol{\Omega}) = I_0\, \boldsymbol{\Omega} \times \boldsymbol{\Omega} = \mathbf{0}$ a rovnice se redukuje na:

$$
I_0\, \dot{\boldsymbol{\Omega}} \;=\; \boldsymbol{\tau}^{\mathrm{body}}.
$$

Toto je rozpojené napříč složkami a identické s druhým Newtonovým zákonem aplikovaným na jediný rotační stupeň volnosti. V bezsilovém případě ($\boldsymbol{\tau} = \mathbf{0}$) je $\boldsymbol{\Omega}$ zachováno v tělesové soustavě a skrze $\dot{\mathbf{R}} = \mathbf{R}\,\widehat{\boldsymbol{\Omega}}$ se inerciální orientace rovnoměrně otáčí kolem konstantního vektoru úhlové rychlosti. Žádná precese; žádná nutace; těleso je rotačně „nudné“. V případě gravitačního gradientu moment síly z gravitačního gradientu mizí (kap. 03 §3.8.2: pro $\mathbf{I}_A = I_0 \mathbb{I}_3$, $\hat{\boldsymbol{\rho}}_{\mathrm{body},A} \times (\mathbf{I}_A \hat{\boldsymbol{\rho}}_{\mathrm{body},A}) = I_0\, \hat{\boldsymbol{\rho}}_{\mathrm{body},A} \times \hat{\boldsymbol{\rho}}_{\mathrm{body},A} = \mathbf{0}$) a dynamika se rozpojí na čistou Keplerovu orbitu + rovnoměrnou inerciální rotaci. Toto je testovací případ §2.4 / kap. 06 limita §6.3 (s analytickou redukcí dále redukující na §6.2 Kepler, je-li $B$ rovněž sférické nebo bodová hmota).

### §5.3.2 Axisymetrický případ — volný symetrický setrvačník

Pro axisymetrické těleso $I_1 = I_2 \equiv I_\perp$, $I_3 \equiv I_\parallel$ v bezsilovém případě se zarámované §5.2.1 rovnice stávají:

$$
I_\perp \dot\Omega_x = (I_\perp - I_\parallel)\,\Omega_y \Omega_z, \qquad
I_\perp \dot\Omega_y = -(I_\perp - I_\parallel)\,\Omega_x \Omega_z, \qquad
I_\parallel \dot\Omega_z = 0.
$$

(Rovnice $\dot\Omega_z$ se zjednodušila, neboť $(I_1 - I_2) = 0$ zabíjí člen vektorového součinu.) Třetí rovnice říká $\Omega_z = \mathrm{const}$ — tělesový spin kolem osy symetrie je zachován. Dosazením zpět do prvních dvou:

$$
\dot\Omega_x = \lambda\, \Omega_y, \qquad \dot\Omega_y = -\lambda\, \Omega_x, \qquad \lambda \equiv \frac{I_\parallel - I_\perp}{I_\perp}\, \Omega_z = \mathrm{const}.
$$

Toto je **Eulerova precese** příčných složek úhlové rychlosti: $\Omega_x, \Omega_y$ oscilují harmonicky s úhlovou frekvencí $\lambda$. Uzavřené řešení je:

$$
\Omega_x(t) = A \cos(\lambda t + \phi_0), \qquad \Omega_y(t) = A \sin(\lambda t + \phi_0), \qquad \Omega_z(t) = \mathrm{const},
$$

s $A$ a $\phi_0$ nastavenými počátečními podmínkami. Vektor tělesové úhlové rychlosti $\boldsymbol{\Omega}(t)$ vykresluje **tělesový kužel** polovičního úhlu $\arctan(A/\Omega_z)$ kolem osy symetrie $\hat{\mathbf{e}}_3$ konstantní rychlostí $\lambda$. Toto je analytické řešení, jež testbed fáze 1 D8 `backends/fortran/tests/test_symmetric_top.f90` validuje numericky, s Podolský 2018 §8.1 Rovn. (8.1)–(8.2) jako přímým CS-jazykovým analytickým orákulem.

### §5.3.3 Asymetrický případ — nestabilita Džanibekova

Pro asymetrické těleso s $I_1 < I_2 < I_3$ v bezsilovém případě mají zarámované §5.2.1 rovnice všechny tři členy vektorového součinu nenulové. Existují tři rotační rovnováhy (rotace kolem každé hlavní osy), se stabilitou klasifikovanou takto:

| Rovnováha | Linearizovaná vlastní hodnota $\sigma^2$ | Klasifikace |
|---|---|---|
| Kolem $\hat{\mathbf{e}}_1$ ($I_1$, nejmenší moment) | $-(\Omega_1^0)^2 (I_2 - I_1)(I_3 - I_1)/(I_2 I_3) < 0$ | Stabilní (oscilační) |
| Kolem $\hat{\mathbf{e}}_2$ ($I_2$, prostřední moment) | $+(\Omega_2^0)^2 (I_2 - I_1)(I_3 - I_2)/(I_1 I_3) > 0$ | **Nestabilní** (exponenciálně narůstající) |
| Kolem $\hat{\mathbf{e}}_3$ ($I_3$, největší moment) | $-(\Omega_3^0)^2 (I_3 - I_1)(I_3 - I_2)/(I_1 I_2) < 0$ | Stabilní (oscilační) |

Nestabilní případ prostřední osy je **efekt tenisové rakety Džanibekova**: malé perturbace spinu srovnaného s prostřední osou narůstají exponenciálně, dokud amplituda nesaturuje, spin se přeskočí přibližně o 180°, pak se re-stabilizuje poblíž opačné konfigurace srovnané s prostřední osou; periodicky poté těleso podstupuje charakteristické přeskoky. Pohyb je omezený (zachování energie + momentu hybnosti jej nutí na uzavřenou trajektorii na kouli tělesové úhlové rychlosti), avšak osa spinu se přerušovaně přeskakuje. Perioda přeskoku závisí na amplitudě počáteční perturbace: větší perturbace → rychlejší přeskoky. Úplné analytické řešení zahrnuje Jacobiho eliptické funkce (Goldstein 2002 §5.6); linearizovaná vlastní hodnota nestability výše dává čas e-násobku raného růstu.

Toto je testbed fáze 1 D9 případ 2 `backends/fortran/tests/test_instability_dzhanibekov.f90`. Každý křížový test fáze 1 chytá neshody znaménka v zarámovaném §5.2.1 tvaru validací ukotvenou ve fyzikální pravdě (přeskok *musí* nastat, perioda *musí* odpovídat predikci eliptické funkce); překlep ve členu vektorového součinu by vystoupil jako přeskok ve špatném směru nebo jako žádný přeskok vůbec.

## §5.4 Dosazení momentu síly z gravitačního gradientu do tvaru hlavních os

Vlož tělesový moment síly z gravitačního gradientu ze zarámovaného výsledku kapitoly 03 §3.7,

$$
\boldsymbol{\tau}_X^{\mathrm{body}} \;=\; \frac{3 G m_{\bar X}}{|\boldsymbol{\rho}|^3}\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, X} \times \big(\mathbf{I}_X \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, X}\big),
$$

kde $\bar X$ je *druhé* těleso, do zarámovaných §5.2.1 Eulerových rovnic hlavních os. V soustavě hlavních os nechť $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, X} = (n_1, n_2, n_3)$ (jednotkový směr mezitělesové vzdálenosti vyjádřený v soustavě hlavních os tělesa $X$, takže $n_1^2 + n_2^2 + n_3^2 = 1$). Pak $\mathbf{I}_X \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},X} = (I_1 n_1, I_2 n_2, I_3 n_3)$ a vektorový součin je:

$$
\hat{\boldsymbol{\rho}}_{\mathrm{body},X} \times (\mathbf{I}_X \hat{\boldsymbol{\rho}}_{\mathrm{body},X})
\;=\;
\begin{vmatrix} \hat{\mathbf{e}}_1 & \hat{\mathbf{e}}_2 & \hat{\mathbf{e}}_3 \\ n_1 & n_2 & n_3 \\ I_1 n_1 & I_2 n_2 & I_3 n_3 \end{vmatrix}
\;=\;
\begin{pmatrix}
(I_3 - I_2)\, n_2 n_3 \\
(I_1 - I_3)\, n_3 n_1 \\
(I_2 - I_1)\, n_1 n_2
\end{pmatrix}.
$$

Moment síly z gravitačního gradientu ve složkách hlavních os je tedy:

$$
\boxed{
\boldsymbol{\tau}_X^{\mathrm{body}}
\;=\;
\frac{3 G m_{\bar X}}{|\boldsymbol{\rho}|^3}\,
\begin{pmatrix}
(I_3^X - I_2^X)\, n_2 n_3 \\
(I_1^X - I_3^X)\, n_3 n_1 \\
(I_2^X - I_1^X)\, n_1 n_2
\end{pmatrix}.
}
$$

Pravá strana nese **tytéž předfaktory typu $(I_i - I_j)$** jako Eulerovy smíšené členy na levé straně zarámovaných §5.2.1 rovnic. To není náhoda: obě pocházejí z antisymetrické / komutátorové struktury $\mathfrak{so}(3)$ působící na symetrické $\mathbf{I}$ a obě identicky mizí pro sférické těleso (všechna $I_i - I_j = 0$).

**Kombinované uzavřené pohybové rovnice hlavních os** (dosazením zarámovaného momentu síly z gravitačního gradientu do zarámovaného §5.2.1 tvaru):

$$
\begin{aligned}
I_1 \dot\Omega_x - (I_2 - I_3)\,\Omega_y\Omega_z \;&=\; (3Gm_{\bar X}/|\boldsymbol{\rho}|^3)\,(I_3 - I_2)\, n_2 n_3, \\
I_2 \dot\Omega_y - (I_3 - I_1)\,\Omega_z\Omega_x \;&=\; (3Gm_{\bar X}/|\boldsymbol{\rho}|^3)\,(I_1 - I_3)\, n_3 n_1, \\
I_3 \dot\Omega_z - (I_1 - I_2)\,\Omega_x\Omega_y \;&=\; (3Gm_{\bar X}/|\boldsymbol{\rho}|^3)\,(I_2 - I_1)\, n_1 n_2.
\end{aligned}
$$

Toto je **uzavřená sada rovnic hlavních os**, již větev Eulerových rovnic `dynamics.py` inženýrské úrovně (tabulka kap. 03 §3.9 → rotační řádky `dynamics.py::state_derivative`) propaguje a již fortranový lege-artis backend `backends/fortran/src/jwst_l2_dynamics.f90` implementuje ve fázi 1.

### §5.4.1 Specializace na librační rovnováhu (řeší OQ-5.2 skrze křížovou kontrolu OQ-FORT-1)

Specializuj na axisymetrické těleso $A$ ($I_1 = I_2 = I_\perp$, $I_3 = I_\parallel$, s hlavními osami srovnanými s tělesovou soustavou) v **radiálně směřující rovnováze**: osa symetrie tělesa srovnaná s $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$, takže $\hat{\boldsymbol{\rho}}_{\mathrm{body},A} = (0, 0, 1)$ → $(n_1, n_2, n_3) = (0, 0, 1)$ v rovnováze. Dosazením do zarámovaného §5.4 momentu síly z gravitačního gradientu:

$$
\boldsymbol{\tau}_A^{\mathrm{body}}\big|_{\mathrm{eq}}
\;=\;
\frac{3Gm_B}{|\boldsymbol{\rho}|^3}\,
\big((I_3 - I_2)\cdot 0 \cdot 1, \; (I_1 - I_3)\cdot 1 \cdot 0, \; (I_2 - I_1)\cdot 0 \cdot 0\big)
\;=\;
\mathbf{0}.
$$

Rovnováha je pravým kritickým bodem: nulový čistý moment síly. **Nyní perturbuj orientaci o malý úhel $\theta$ kolem, řekněme, tělesové osy $\hat{\mathbf{e}}_1$.** Přijmi konvenci kap. 03 §3.7 pro infinitezimální tělesovou rotaci: $\mathbf{R}_q^A \to \mathbf{R}_q^A(\mathbb{I} + \widehat{\delta\boldsymbol{\theta}_A})$ s $\delta\boldsymbol{\theta}_A = (\theta, 0, 0)$. Pak identitou kap. 03 §3.7 $\delta\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = -\delta\boldsymbol{\theta}_A \times \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$:

$$
\delta\hat{\boldsymbol{\rho}}_{\mathrm{body},A}
\;=\;
-(\theta, 0, 0) \times (0, 0, 1)
\;=\;
-\theta\,(\hat{\mathbf{e}}_1 \times \hat{\mathbf{e}}_3)
\;=\;
-\theta\,(-\hat{\mathbf{e}}_2)
\;=\;
+\theta\,\hat{\mathbf{e}}_2.
$$

Perturbovaný tělesový směr vzdálenosti je tedy $\hat{\boldsymbol{\rho}}_{\mathrm{body},A}(\theta) \approx (0, +\theta, 1)$, dávající $(n_1, n_2, n_3) = (0, +\theta, 1)$, a perturbovaný moment síly je:

$$
\boldsymbol{\tau}_A^{\mathrm{body}}(\theta)
\;\approx\;
\frac{3Gm_B}{|\boldsymbol{\rho}|^3}\,
\big((I_3 - I_2)\cdot (+\theta) \cdot 1, \; 0, \; 0\big)
\;=\;
\frac{3Gm_B(I_\parallel - I_\perp)}{|\boldsymbol{\rho}|^3}\,(\theta, 0, 0),
$$

s užitím $I_2 = I_\perp, I_3 = I_\parallel$ pro axisymetrické těleso. Librační rovnice pro perturbaci malého úhlu je pak $I_\perp\, \ddot\theta = \tau_x^{\mathrm{body}}(\theta)$:

$$
\boxed{
I_\perp\, \ddot\theta \;=\; \frac{3Gm_B(I_\parallel - I_\perp)}{|\boldsymbol{\rho}|^3}\, \theta.
}
$$

Toto je lineární ODR $\ddot\theta = k\, \theta$ s $k = 3Gm_B(I_\parallel - I_\perp)/(I_\perp\,|\boldsymbol{\rho}|^3)$. Znaménko $k$ určuje stabilitu:

- **Protáhlé** ($I_\parallel < I_\perp$, $k < 0$): $\ddot\theta = -|k|\theta$ → jednoduchý harmonický oscilátor → **stabilní librace** s frekvencí $\omega_{\mathrm{lib}}^2 = |k| = 3Gm_B(I_\perp - I_\parallel)/(I_\perp\,|\boldsymbol{\rho}|^3) > 0$. Postoj s osou symetrie podél radiály je gravitačně-gradientně stabilizovaným postojem (těleso „chce“ mířit svou dlouhou osou na primár).
- **Zploštělé** ($I_\parallel > I_\perp$, $k > 0$): $\ddot\theta = +k\theta$ → exponenciálně narůstající → **nestabilní** s časem e-násobku $\tau_{\mathrm{inst}} = 1/\sqrt{k}$. Stabilní postoj má osu symetrie *kolmou* k $\hat{\boldsymbol{\rho}}$ (ploché těleso podobné sluneční cloně se srovnává čelem).

**Srovnání s orientační Jacobiho maticí kap. 03 §3.7.2.** Kap. 03 §3.7.2 zarámovaná (po opravě znaménka OQ-FORT-1 2026-05-24 odp.) dává matici tuhosti librační rovnováhy:

$$
\mathbf{J}^{\mathrm{orient}}_{A}\big|_{\mathrm{eq}}
\;=\; \alpha\, (I_\parallel - I_\perp)\, \mathrm{diag}(1, 1, 0)
\;=\; \frac{3Gm_B(I_\parallel - I_\perp)}{|\boldsymbol{\rho}|^3}\,\mathrm{diag}(1,1,0),
$$

s vlastní hodnotou $J_{11} = 3Gm_B(I_\parallel - I_\perp)/|\boldsymbol{\rho}|^3$. Vynásobením linearizovaného vztahu $\delta\boldsymbol{\tau}_A^{\mathrm{body}} = \mathbf{J}^{\mathrm{orient}}_A\, \delta\boldsymbol{\theta}_A$ do librační rovnice $I_\perp\ddot\theta = \delta\tau_x = J_{11}\theta$:

$$
I_\perp\,\ddot\theta \;=\; J_{11}\, \theta \;=\; \frac{3Gm_B(I_\parallel - I_\perp)}{|\boldsymbol{\rho}|^3}\,\theta,
$$

což je **identické** se zarámovaným výsledkem tohoto oddílu. Obě cesty — (a) dosazení zarámovaného §5.4 momentu síly ve tvaru hlavních os do zarámované Eulerovy rovnice kap. 05 §5.2.1 v librační rovnováze (tato §5.4.1) a (b) aplikace orientační Jacobiho matice kap. 03 §3.7.2 v rovnováze přímo — produkují tutéž linearizovanou librační rovnici s toutéž klasifikací stability (protáhlé stabilní / zploštělé nestabilní). **(Řeší OQ-5.2.)**

Toto je **disciplína dosazení pro kontrolu zdravého rozumu zarámovaného vzorce** (dle KB-OPUS-DOC-INTRA-INCONSISTENCY): §5.4 zarámovaný moment síly z gravitačního gradientu ve tvaru hlavních os je ověřen vůči nezávisle odvozené kap. 03 §3.7.2 zarámované orientační Jacobiho matici dosazením specifického vypracovaného případu (axisymetrické těleso, rovnováha na ose, příčná perturbace) a potvrzením shody. Testbed fáze 1 D5 `backends/fortran/tests/test_orientation_jacobian.f90` (T2a librační-rovnovážné vlastní hodnoty) je regresní pojistkou mezi implementacemi při toleranci ULP × √N — táž data, dvě odvozovací cesty, tři implementace (kanonická úroveň kap. 03, kanonická úroveň kap. 05, inženýrská úroveň Fortran D5) vše ve shodě.

### §5.4.2 Struktura tříd symetrie dosazených rovnic

Čtením tabulky tříd symetrie z §5.3 ve spojení se zarámovaným §5.4 momentem síly z gravitačního gradientu:

| Třída symetrie tělesa | Počet nenulových složek momentu síly hlavních os | Dynamický obsah |
|---|---|---|
| Sférická ($I_1 = I_2 = I_3$) | 0 (všechna $(I_i - I_j) = 0$) | Moment síly z gravitačního gradientu identicky mizí; rotace se rozpojí od orbity; testbed §2.4 / kap. 06 limita §6.3. |
| Axisymetrická ($I_1 = I_2 \neq I_3$) | 2 genericky ($\tau_z^{\mathrm{body}}$ mizí, když $I_1 = I_2$, ponechávajíc $\tau_x, \tau_y$ jako librační-hnací složky) | Dvourozměrná librační rovina kolmá k ose symetrie; jeden cyklický směr podél osy. Testbed §2.5 / kap. 06 §6.4 (librace malého úhlu / Lagrangeův setrvačník). |
| Asymetrická (všechna $I_i$ odlišná) | 3 genericky | Tři nezávislé librační / nestabilní módy; moment síly z gravitačního gradientu se váže na všechny tři orientační stupně volnosti. Vícemódová vazba; ne-uzavřená kromě zvláštních případů. Testbed §2.5 + §3.5 / kap. 06 §6.4-rozšířený. |

Struktura je **zrcadlově symetrická** s tabulkou tříd symetrie z §5.3 pro bezsilový případ: táž podmínka anizotropie $(I_i - I_j) \neq 0$, jež činí setrvačný smíšený člen aktivním v bezsilové Eulerově rovnici, činí složku momentu síly z gravitačního gradientu aktivní v dosazené rovnici. To není náhoda: oba členy pocházejí z antisymetrické tenzorové struktury $\boldsymbol{\Omega} \times (\mathbf{I}\boldsymbol{\Omega})$ vs $\hat{\boldsymbol{\rho}}_{\mathrm{body}} \times (\mathbf{I}\hat{\boldsymbol{\rho}}_{\mathrm{body}})$ a mají paralelní tabulky hodností.

## §5.5 Integrály pohybu a strukturální zákony zachování

Tělesové Eulerovy rovnice mají v bezsilovém případě dobře známé integrály pohybu:

- **Kinetická energie v tělesové soustavě** $T = \tfrac{1}{2}(I_1 \Omega_x^2 + I_2 \Omega_y^2 + I_3 \Omega_z^2)$ je zachována.
- **Velikost tělesového momentu hybnosti** $|\mathbf{L}|^2 = (I_1 \Omega_x)^2 + (I_2 \Omega_y)^2 + (I_3 \Omega_z)^2$ je zachována.

**Ověření.** Derivuj $T$ podle času:

$$
\frac{dT}{dt} \;=\; I_1 \Omega_x \dot\Omega_x + I_2 \Omega_y \dot\Omega_y + I_3 \Omega_z \dot\Omega_z.
$$

Dosazením bezsilových zarámovaných §5.2.1 Eulerových rovnic ($I_i \dot\Omega_i = (I_j - I_k)\Omega_j\Omega_k$ pro příslušné cyklické indexy):

$$
\frac{dT}{dt}
\;=\;
\Omega_x \cdot (I_2 - I_3)\Omega_y\Omega_z + \Omega_y \cdot (I_3 - I_1)\Omega_z\Omega_x + \Omega_z \cdot (I_1 - I_2)\Omega_x\Omega_y
\;=\;
\Omega_x \Omega_y \Omega_z \cdot \big[(I_2 - I_3) + (I_3 - I_1) + (I_1 - I_2)\big]
\;=\;
0.
$$

Tři koeficienty se identicky sčítají na nulu, potvrzujíc $T = \mathrm{const}$. Analogický výpočet pro $|\mathbf{L}|^2 = (I_1\Omega_x)^2 + (I_2\Omega_y)^2 + (I_3\Omega_z)^2$ postupuje:

$$
\frac{d|\mathbf{L}|^2}{dt} = 2\, I_1\Omega_x \cdot I_1 \dot\Omega_x + \ldots = 2\, I_1\Omega_x \cdot (I_2 - I_3)\Omega_y\Omega_z + \ldots = 2\Omega_x\Omega_y\Omega_z\, [I_1(I_2 - I_3) + I_2(I_3 - I_1) + I_3(I_1 - I_2)] = 0,
$$

opět algebraickou identitou $I_1(I_2 - I_3) + I_2(I_3 - I_1) + I_3(I_1 - I_2) = 0$ (rozviň a vyruš).

### §5.5.1 Geometrická interpretace

Bezsilová trajektorie v tělesovém $\boldsymbol{\Omega}$-prostoru je **průnikem** dvou elipsoidů:

- **Elipsoidu kinetické energie** $\tfrac{1}{2}(I_1 \Omega_x^2 + I_2 \Omega_y^2 + I_3 \Omega_z^2) = T_0$ (pevná energetická plocha).
- **Elipsoidu velikosti momentu hybnosti** $(I_1\Omega_x)^2 + (I_2\Omega_y)^2 + (I_3\Omega_z)^2 = L_0^2$ (pevná plocha čtverce momentu hybnosti).

Oba jsou hladké elipsoidy v $\boldsymbol{\Omega}$-prostoru. Jejich průnikem je 1rozměrná křivka (genericky; degenerované případy ve zvláštních pevných bodech tříd symetrie). Trajektorie $\boldsymbol{\Omega}(t)$ vykresluje tuto křivku; je uzavřená (neboť obě plochy jsou omezené) a periodická kromě nestabilních pevných bodů prostřední osy z §5.3.3, kde křivka prochází sedlem. Uzavřené řešení zahrnuje Jacobiho eliptické funkce; analytický detail je v Landau-Lifšic, sv. I §37 a Goldstein 2002 §5.6. **(Řeší OQ-5.3 částečně — bezsilové konstanty jsou stanoveny; externě-momentovaný případ je nyní pojednán.)**

### §5.5.2 S vnějším momentem síly — zachování celkové energie

Když těleso zakouší vnější moment síly (náš případ gravitačního gradientu), kinetická energie $T$ **není** zachována sama. Věta o práci a energii dává:

$$
\frac{dT}{dt} \;=\; \boldsymbol{\Omega} \cdot \boldsymbol{\tau}^{\mathrm{body}},
$$

což je výkon dodaný rotačnímu stupni volnosti vnějším momentem síly. (Odvození: derivuj $T = \tfrac{1}{2}\boldsymbol{\Omega}^T \mathbf{I} \boldsymbol{\Omega}$; dosaď tělesovou Eulerovu rovnici $\mathbf{I}\dot{\boldsymbol{\Omega}} = -\boldsymbol{\Omega}\times(\mathbf{I}\boldsymbol{\Omega}) + \boldsymbol{\tau}^{\mathrm{body}}$; užij cyklickou identitu smíšeného součinu $\boldsymbol{\Omega}\cdot[\boldsymbol{\Omega}\times(\mathbf{I}\boldsymbol{\Omega})] = 0$.)

Moment síly je tělesovou projekcí $-\partial V/\partial \boldsymbol{\theta}$; odpovídající **časová derivace potenciální energie** je:

$$
\frac{dV_{\mathrm{tidal},X}}{dt}
\;=\;
\frac{\partial V_{\mathrm{tidal},X}}{\partial q_X} \cdot \dot q_X
\;+\; \frac{\partial V_{\mathrm{tidal},X}}{\partial \boldsymbol{\rho}} \cdot \dot{\boldsymbol{\rho}}
\;=\;
-\boldsymbol{\tau}_X^{\mathrm{body}} \cdot \boldsymbol{\Omega}_X
\;+\; \nabla_{\boldsymbol{\rho}} V_{\mathrm{tidal},X}\cdot \dot{\boldsymbol{\rho}},
$$

s prvním členem jako rychlostí, kterou moment síly koná zápornou práci na tělese (těleso získává $T$ rychlostí, kterou orientační potenciál ztrácí $V$), a druhým jako rychlostí, kterou zpětná síla (kap. 04 §4.6.2 zarámovaná) koná práci na relativní orbitě. Sečtením:

$$
\frac{d(T + V)}{dt} \;=\; 0 \quad \text{(zachování celkové energie)}.
$$

Úplným zákonem zachování je tedy **celková mechanická energie** $E = T + V$, s $V$ daným kap. 03 §3.4.3 zarámovaným useknutým vzájemným potenciálem. Newtonovsko-třetí-zákon rovnováha kap. 04 §4.6.4 plus tento rotační vztah práce-energie dávají:

- $\boldsymbol{P}_{\mathrm{CoM}}$ zachováno (translační invariance — kap. 02 §2.6.1 + Noether).
- $\boldsymbol{L}_{\mathrm{total}}$ zachováno (rotační invariance — kap. 02 + kap. 04 §4.6.4).
- $E = T + V$ zachováno (invariance vůči časové translaci + účetnictví práce-energie / potenciální-energie výše).

Toto jsou **tři strukturální zákony zachování**, jež musí jakákoli poctivá integrace zarámované §4.7 sady rovnic respektovat. Diagnostiky inženýrské úrovně z první verze hlásí $|dE/E_0| = 6{,}2 \times 10^{-12}$ a $|d\mathbf{L}/L_0| = 2{,}1 \times 10^{-10}$ přes 600 sekund při $dt = 0{,}05$ s — malý, avšak nenulový úlet je kumulativní integrační chybou RK4 (nesymplektickou; škálování $O(dt^4 \cdot T)$), nikoli porušením strukturálních zákonů zachování. **(Řeší OQ-5.3 plně — bezsilové konstanty $(T, |\mathbf{L}|^2)$ i momentované zachování celkové energie $E = T + V$ obě stanoveny s explicitním účetnictvím práce-energie.)**

## §5.6 Křížové odkazy na validační testovací případy (ukazatele na kapitolu 06)

Každá specializace třídy symetrie z §5.3 má řešení analytického limitu, jež se stává validačním testovacím případem v kapitole 06:

- **Sférické těleso, bezsilové** → triviální validace: $\boldsymbol{\Omega} = \mathrm{const}$, orientace rotuje rovnoměrně. Kap. 06 §6.3-sférická-limita / `_specs/...` §3.4 rozšířený případ počátečních podmínek nulové rotace.
- **Axisymetrické těleso, bezsilové** → **volný symetrický setrvačník** s regulární precesí; testbed validuje tělesové řešení vůči analytickému vzorci $\Omega_x(t) = A \cos(\lambda t + \phi_0)$, $\Omega_y(t) = A \sin(\lambda t + \phi_0)$, $\Omega_z = \mathrm{const}$. Toto je test fáze 1 D8 `backends/fortran/tests/test_symmetric_top.f90`, ukotvený k Podolský §8.1 Rovn. (8.1)–(8.2) jako CS-jazykové analytické orákulum. Kap. 06 §6.3.
- **Asymetrické těleso, bezsilové** → **efekt tenisové rakety Džanibekova**; testbed validuje predikované přeskoky o konečné amplitudě. Test fáze 1 D9 případ 2 `backends/fortran/tests/test_instability_dzhanibekov.f90`. Kap. 06 §6.5.
- **Axisymetrické těleso, moment síly z gravitačního gradientu (malá librace)** → linearizovaná librační tuhost odpovídá výsledku kapitoly 03 §3.7.2 $\alpha(I_\parallel - I_\perp)$; trajektorie plné amplitudy odpovídá analýze efektivního potenciálu **Lagrangeova setrvačníku**. Testbed fáze 2; ukotvený k Podolský §8.3 Rovn. (8.10)–(8.15) jako CS-jazykové analytické orákulum. Kap. 06 §6.4.
- **Asymetrické těleso, moment síly z gravitačního gradientu** → vícemódová librace se stabilitou klasifikovanou plnou vlastní strukturou Jacobiho matice (kapitola 03 §3.7.2). Testbed fáze 2; žádné uzavřené analytické řešení, validuje skrze bitovou identitu mezi implementacemi (doktrína §4.4 C2). Kap. 06 §6.5-rozšířený.

Kapitola 06 přivádí každý testbed ke kvantitativní bráně vyhověl/nevyhověl.

## §5.7 Co je stanoveno na konci této kapitoly

Čtenář, jenž prošel §5.1–§5.6, má:

- Klasické Eulerovy dynamické rovnice z roku 1758 (zarámované §5.2.1) odvozené rozvojem složka po složce souřadnicově nezávislého tělesového Lagrangiánu odvozeného Euler-Poincaré kapitoly 04 §4.4 — potvrzujíc strukturální specializaci. **OQ-5.1 vyřešeno.**
- Kontrolu zdravého rozumu konvence znaménka (§5.2.2) skrze vlastní hodnotu nestability Džanibekova, konzistentní s učebnicovou klasifikací teorému tenisové rakety.
- Explicitní **specializace tříd symetrie** (§5.3): sférickou (žádná vnitřní rotační dynamika), axisymetrickou (volný symetrický setrvačník s uzavřeným precesním řešením při polovičním úhlu tělesového kužele $\arctan(A/\Omega_z)$ a precesní rychlosti $\lambda = (I_\parallel - I_\perp)\Omega_z/I_\perp$) a asymetrickou (nestabilita Džanibekova kolem prostřední hlavní osy s explicitními vzorci vlastních hodnot ve třech rovnováhách hlavních os).
- **Moment síly z gravitačního gradientu dosazený explicitně** do tvaru hlavních os (§5.4 zarámovaný); pozorována strukturální zrcadlová symetrie mezi setrvačně-anizotropní strukturou bezsilových Eulerových smíšených členů a strukturou složek momentu síly z gravitačního gradientu.
- **Specializaci na librační rovnováhu** (§5.4.1) reprodukující zarámovanou orientační Jacobiho matici kap. 03 §3.7.2 s vlastní hodnotou $J_{11} = 3Gm_B(I_\parallel - I_\perp)/|\boldsymbol{\rho}|^3$ nezávislou odvozovací cestou. **OQ-5.2 vyřešeno** skrze disciplínu dosazení pro kontrolu zdravého rozumu zarámovaného vzorce; librační analýza z dosazených rovnic této kapitoly souhlasí s analýzou orientační Jacobiho matice kap. 03 §3.7.2. Regresní pojistka mezi implementacemi u testbedu fáze 1 D5.
- **Strukturu integrálů pohybu** (§5.5.1 + §5.5.2): v bezsilovém případě jsou tělesová kinetická energie a čtverec momentu hybnosti zachovány identitou $I_1(I_2-I_3) + I_2(I_3-I_1) + I_3(I_1-I_2) = 0$; geometrická interpretace jako průnik energetického elipsoidu a elipsoidu velikosti momentu hybnosti v tělesovém $\boldsymbol{\Omega}$-prostoru. V externě-momentovaném případě je celková mechanická energie $E = T + V$ zachována, s explicitním účetnictvím práce-energie / potenciální-energie. **OQ-5.3 vyřešeno.**
- **Mapu ukazatelů** z každé třídy symetrie na odpovídající validační testovací případ analytického limitu kapitoly 06 (§5.6).

**Co dosud stanoveno NENÍ:**

- Úplná validace vůči testovacím případům analytických limitů při kvantitativních branách vyhověl/nevyhověl. **Kapitola 06.**
- Symplektický integrátor na Lieově grupě, jenž by nahradil RK4 pro dlouhohorizontovou přesnou integraci těchto rovnic. **Rozsah inženýrské úrovně v0.2** (kap. 04 §4.9 OQ-4.6).
- Rozšíření L2 / omezené soustavy tří těles přidávající centrifugální a Coriolisovy pseudo-síly k výrazu momentu síly z gravitačního gradientu. **Kanonická úroveň v0.2.**
- Plná nelineární analýza stability v libračních rovnováhách energeticko-momentovou metodou. **Rozsah kanonické úrovně v0.2** (kap. 04 §4.9 OQ-4.7).

---

## §5.8 Otevřené otázky v této kapitole

### Vyřešeno v tomto kole (kanonická úroveň, 3. kolo, 2026-05-25)

- **OQ-5.1** ***Vyřešeno*** — v §5.2.1: explicitní rozvoj složka po složce $\boldsymbol{\Omega}_X \times (\mathbf{I}_X \boldsymbol{\Omega}_X)$ v soustavě hlavních os produkuje zarámované skalární Eulerovy rovnice $I_1 \dot\Omega_x - (I_2 - I_3)\Omega_y\Omega_z = \tau_x^{\mathrm{body}}$ (a cyklicky). Obnoveno z kap. 04 §4.4 souřadnicově nezávislého výsledku strukturální specializací.
- **OQ-5.2** ***Vyřešeno*** — v §5.4.1: specializace na librační rovnováhu produkuje librační frekvenci $\omega^2 = 3Gm_B(I_\parallel - I_\perp)/(I_\perp |\boldsymbol{\rho}|^3)$ přímou linearizací zarámovaných §5.4 dosazených Eulerových rovnic. Velikost tuhosti odpovídá zarámované orientační Jacobiho matici kap. 03 §3.7.2 s vlastní hodnotou $J_{11} = \alpha(I_\parallel - I_\perp)$ přesně. Konvence znaménka (protáhlé stabilní / zploštělé destabilizující) konzistentní s tabulkou stability kap. 03 §3.7.2 po opravě znaménka OQ-FORT-1. Disciplína dosazení pro kontrolu zdravého rozumu zarámovaného vzorce potvrzena shodou mezi nezávislými odvozeními kap. 03 a kap. 05.
- **OQ-5.3** ***Vyřešeno*** — v §5.5.1 + §5.5.2: bezsilové konstanty $T$ a $|\mathbf{L}|^2$ zachovány explicitní algebraickou identitou (součet cyklů $(I_i - I_j)$ roven nule); momentovaný-případ celková energie $E = T + V$ zachována účetnictvím práce-energie / potenciální-energie, s explicitním odvozením v §5.5.2. Diagnostiky zachování fáze 1 $|dE/E_0|, |d\mathbf{L}/L_0|$ jsou regresní pojistkou inženýrské úrovně.
- **OQ-5.4** ***Vyřešeno*** — v §5.6 + křížová kontrola vůči tabulce kap. 04 §4.7.1: zarámované §5.4 dosazené rovnice hlavních os jsou řádek po řádku konzistentní s rotačními řádky `dynamics.py::state_derivative` (162–163) + propagací Eulerových rovnic fortranového `backends/fortran/src/jwst_l2_dynamics.f90` z fáze 1. Mezi­úrovňová shoda je přesná na rotační úrovni; jedinou mezerou je zpětná síla na orbitu (§5.4 dosazené rovnice popisují pouze rotační dynamiku; zpětná síla je územím kap. 04 §4.6.2 + §4.7.1 / §4.9 OQ-4.5, nikoli kap. 05).

### Vyřešeno v rozšíření tohoto kola (kanonická úroveň, 3. kolo+, 2026-05-25)

- **OQ-5.5** ***Vyřešeno*** — v kap. 06 §6.5.3: uzavřené Jacobiho-eliptické řešení bezsilového problému asymetrického setrvačníku zpracováno. Zarámovaná parametrizace $(L_x, L_y, L_z) = (\alpha\sqrt{\cdot}\,\mathrm{cn}(\Omega^*t, k), \alpha\,\mathrm{sn}(\Omega^*t, k), \alpha\sqrt{\cdot}\,\mathrm{dn}(\Omega^*t, k))$ s charakteristickou frekvencí $\Omega^*$ a eliptickým modulem $k$ odvozenými z $T$ a $|\mathbf{L}|^2$. Zarámovaná perioda $T_{\mathrm{LL}} = 4K(k)/\Omega^*$ → perioda přeskoku Džanibekova jako falzifikovatelná funkce amplitudy počáteční perturbace. Hladká interpolace mezi §6.3 volným symetrickým setrvačníkem ($k \to 0$) a nestabilní rovnováhou prostřední osy ($k \to 1$). Anchor: Landau-Lifšic, sv. I §37, Goldstein 2002 §5.6.
- **OQ-5.6** ***Vyřešeno*** — v kap. 06 §6.4.3: uzavřené řešení Lagrangeova setrvačníku + librace gravitačního gradientu metodou efektivního potenciálu $V_{\mathrm{ef}}(\vartheta)$. Zarámovaný efektivní potenciál $V_{\mathrm{ef}}(\vartheta) = (L_\varphi - L_\psi\cos\vartheta)^2/(2I_\perp\sin^2\vartheta) - 3Gm_B(I_\perp - I_\parallel)\cos^2\vartheta/(2|\boldsymbol{\rho}|^3)$. První integrály cyklických souřadnic $L_\psi, L_\varphi$ odvozeny; uzavřená kvadratura pro $\vartheta(t)$ dána; tři precesní režimy (a)/(b)/(c) klasifikovány dle Podolský §8.3 fenomenologie. Librační perioda konečné amplitudy $T_{\mathrm{lib}}(\vartheta_{\mathrm{amp}}) = (2/\pi)K(k)\cdot T_{\mathrm{lib}}^{\mathrm{linear}}$ jako korekce eliptickou-funkcí-amplitudy. Anchor: Podolský 2018 §8.3 Rovn. (8.10)–(8.15) přímo; Goldstein 2002 §5.7.

---

> **Doplňující anglojazyčné odkazy.** Euler 1758 — *Théorie du mouvement des corps solides ou rigides* (Petrohradská akademie) — původní odvození tělesových rovnic nesoucích jeho jméno. Moderní prezentace: Goldstein 2002 §5.4 (učebnicové pojednání s diskusí hlavních os); Landau-Lifšic 1976 §35 (stručné tenzorové odvození); Arnold 1989 §29 (moderní perspektiva Lieových grup, spojující se s Euler-Poincarého zarámováním kapitoly 04); Hughes 2004 §3.3 (pojednání domény postoje kosmické lodě, s explicitními specializacemi tříd symetrie a analýzou librace s hlavní-osou-srovnanou-s-radiálou, jež ukotvuje §5.6).
>
> **Hlavní odkaz.** Podolský 2018 Kapitola 7 §7.2 *Eulerovy dynamické rovnice* je přímým CS-jazykovým analogem centrálního zarámovaného výsledku této kapitoly (§5.2.1): jeho **Rovn. (7.15) je bit-identická s našimi §5.2.1 zarámovanými Eulerovými rovnicemi**, odvozená v jeho prezentaci cestou Newtonovsko-Eulerovou $d\mathbf{L}/dt = \mathbf{M}$ plus transportní rovnicí (6.8). Obě cesty (Newtonovsko-Eulerova v Podolský §7.2; Euler-Poincarého redukce v naší kapitole 04) produkují tytéž rovnice — týž vzor shoda-jako-validace jako v naší §3.6 a §4.5. Pro spojení s validačními testovacími případy (§5.6) je Podolský Kapitola 8 *Aplikace: setrvačníky* přímým CS-jazykovým odkazem: §8.1 volný symetrický setrvačník + §8.2 setrvačník se třením (tlumený) + §8.3 Lagrangeův setrvačník s metodou $V_{\mathrm{ef}}$. Každý test fáze 1 / fáze 2 v našich portech backendu má své CS-jazykové analytické orákulum v Podolský Kap. 8.

---

**Následující kapitola:** [06 — Analytické limity a validační brány](06-analyticke-limity-a-validacni-brany.md)

**Konec 05-telesove-eulerovy-rovnice.md (v0.1, 3. kolo).**
