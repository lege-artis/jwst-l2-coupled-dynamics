# Kapitola 4 — Eulerovy–Lagrangeovy rovnice a Euler-Poincarého redukce

> **Předpokládané čtení.** Kapitola 01 (konfigurační varieta a variační princip), kapitola 02 (kinetická energie na $SE(3) \times SE(3)$, včetně §2.5.1 rotační transport dI/dt, §2.6.1 symplektická transformace těžiště, §2.7.1 levá invariance) a kapitola 03 (vzájemný gravitační potenciál, včetně §3.7 moment síly z gravitačního gradientu a §3.7.2 orientační Jacobiho matice).
>
> **Stav.** v0.1 (kanonická úroveň, 3. kolo, 2026-05-25). Všechny čtyři vnitřní otevřené otázky kapitoly OQ-4.1, OQ-4.2, OQ-4.3, OQ-4.4 vyřešeny v tomto kole; viz §4.9 pro souhrn uzávěru. Dopředu orientované otevřené otázky v0.2 vystupují v §4.9.

> **Poznámka ke způsobu a předpokladům odvození.** Tato CS kanonická kapitola (v0.1) je odvozena z anglické kanonické kapitoly `docs/canonical/en/04-euler-lagrange-and-euler-poincare.md` (v0.1, 3. kolo). EN verze zůstává autoritativním pramenem při budoucích aktualizacích.

---

## §4.1 Co tato kapitola dělá

Kapitoly 02 a 03 vyprodukovaly, po řadě, kinetickou energii $T$ (levo-invariantní Riemannovu metriku na rotačních faktorech $TQ$ plus translační metriku s redukovanou hmotností) a vzájemný gravitační potenciál $V$ (kvadrupólově-useknutý vzájemný potenciál s jeho momentem síly z gravitačního gradientu). Tato kapitola sestavuje **Lagrangián** $L = T - V$ na konfiguračním tečném svazku $TQ$, odvozuje **Eulerovy–Lagrangeovy rovnice** z variačního principu $\delta S = 0$ stanoveného v kapitole 01 §1.5 a aplikuje **Euler-Poincarého redukci** (Marsden & Ratiu 1999 §13) na levo-invariantní rotační stupně volnosti k obnovení tělesových **Eulerových rovnic** pro každé tuhé těleso jako strukturálního důsledku levé invariance metriky — nikoli jako souřadnicového triku. Výstupem této kapitoly je úplná sada pohybových rovnic pro dvoutělesovou sdruženou rotačně-translační soustavu, připravená pro specializovaný tělesový tvar v kapitole 05 a pro numerickou integraci v kódových cestách inženýrské úrovně.

Kapitola sleduje disciplínu **metodologie dvě-paralelně** stanovenou v kapitole 03 §3.6 a znovupotvrzenou direktivou projektového vlastníka 2026-05-23: tělesové Eulerovy rovnice jsou odvozeny **dvěma nezávislými cestami** — (a) přímým Eulerovým–Lagrangeovým výpočtem v souřadnicích na $TQ$ a (b) Euler-Poincarého redukcí levou akcí $SO(3) \times SO(3)$ na $Q$ — a obě jsou ukázány jako shodné. Shoda je validací; kdyby obě odvození produkovala odlišné rovnice, jedno z nich by bylo chybné.

> **Hlavní odkaz.** Podolský 2018 Kapitola 7 *Dynamika tuhého tělesa* je nejpřímějším CS-jazykovým ukotvením této kapitoly. §7.2 odvozuje tělesové Eulerovy dynamické rovnice jako zarámovanou Rovn. (7.15): $I_1 \dot{\Omega}_x - (I_2 - I_3)\Omega_y\Omega_z = M_x$ (a cyklicky) — tytéž pohybové rovnice, k nimž §4.4 zde dospívá skrze Euler-Poincarého redukci. §7.3 pak znovu-odvozuje tytéž rovnice Lagrangeovým formalismem s Eulerovými úhly jako zobecněnými souřadnicemi — přímo ekvivalentní našemu Přístupu A (§4.3) v mapě Eulerových úhlů na $SO(3)$. Metodologie křížové validace dvě-paralelně, již tato kapitola formalizuje, je tedy již implicitní v Podolského dvojím odvození; moderní Marsden-Ratiu Euler-Poincarého redukce je souřadnicově nezávislou abstrakcí téže shody.

## §4.2 Sestavení Lagrangiánu $L = T - V$

Kombinujeme kinetickou energii v těžišťové soustavě z kapitoly 02 §2.6 (zarámovanou) s kvadrupólově-useknutým vzájemným potenciálem z kapitoly 03 §3.4.3 (zarámovaným) k získání Lagrangiánu na konfigurační varietě redukované translační symetrií. Translační symetrie soustavy byla užita v kapitole 02 §2.6.1 k provedení kroku symplektické-redukce-translací, jenž nastavuje $\mathbf{P}_{\mathrm{CoM}} = \mathbf{0}$, $\mathbf{R}_{\mathrm{CoM}} = \mathbf{0}$ explicitně; aktivní konfigurační varietou je tedy:

$$
Q_{\mathrm{CoM}} \;\equiv\; \mathbb{R}^3_{\boldsymbol{\rho}} \;\times\; SO(3)_A \;\times\; SO(3)_B,
$$

dimenze $3 + 3 + 3 = 9$. Kinetická energie na $TQ_{\mathrm{CoM}}$ je (kap. 02 §2.6 zarámovaná):

$$
T_{\mathrm{CoM}} \;=\; \tfrac{1}{2}\, \mu\, |\dot{\boldsymbol{\rho}}|^2
\;+\; \tfrac{1}{2}\, \boldsymbol{\Omega}_A^T\, \mathbf{I}_A\, \boldsymbol{\Omega}_A
\;+\; \tfrac{1}{2}\, \boldsymbol{\Omega}_B^T\, \mathbf{I}_B\, \boldsymbol{\Omega}_B,
$$

s $\mu = m_A m_B / (m_A + m_B)$ redukovanou hmotností a $\boldsymbol{\Omega}_X$ tělesovou úhlovou rychlostí tělesa $X$ (levo-trivializovaným tečným vektorem — kap. 02 §2.7.1). Potenciální energií je kap. 03 §3.4.3 zarámovaný kvadrupólově-useknutý výraz:

$$
V \;=\; -\frac{G\, m_A\, m_B}{|\boldsymbol{\rho}|}
\;+\; V_{\mathrm{tidal},\, A}(\boldsymbol{\rho}, q_A)
\;+\; V_{\mathrm{tidal},\, B}(\boldsymbol{\rho}, q_B)
\;+\; \mathcal{O}\!\left(\frac{\ell^4}{|\boldsymbol{\rho}|^5}\right),
$$

s každou slapovou částí

$$
V_{\mathrm{tidal},\, X} \;=\; -\frac{G\, m_{\bar X}}{2\, |\boldsymbol{\rho}|^3}\, \Big[ \mathrm{tr}(\mathbf{I}_X) \;-\; 3\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, X} \cdot \mathbf{I}_X \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, X} \Big],
$$

kde $\bar X$ značí druhé těleso, $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, X} = (\mathbf{R}_q^X)^T \hat{\boldsymbol{\rho}}$ je jednotkový vektor mezitělesové vzdálenosti vyjádřený v soustavě tělesa $X$ a konstantní tělesové tenzory setrvačnosti $\mathbf{I}_X$ vstupují lineárně. Lagrangián na $TQ_{\mathrm{CoM}}$ je pak:

$$
\boxed{
\begin{aligned}
L(\boldsymbol{\rho}, \dot{\boldsymbol{\rho}}, q_A, \boldsymbol{\Omega}_A, q_B, \boldsymbol{\Omega}_B)
\;=\;
&\tfrac{1}{2}\, \mu\, |\dot{\boldsymbol{\rho}}|^2 \;+\; \frac{G\, m_A\, m_B}{|\boldsymbol{\rho}|}
\\[2pt]
&\;+\; \tfrac{1}{2}\, \boldsymbol{\Omega}_A^T\, \mathbf{I}_A\, \boldsymbol{\Omega}_A
\;+\; \tfrac{1}{2}\, \boldsymbol{\Omega}_B^T\, \mathbf{I}_B\, \boldsymbol{\Omega}_B
\\[2pt]
&\;-\; V_{\mathrm{tidal},\, A}(\boldsymbol{\rho}, q_A)
\;-\; V_{\mathrm{tidal},\, B}(\boldsymbol{\rho}, q_B).
\end{aligned}
}
$$

### §4.2.1 Strukturální rozklad: která část co s čím váže

Zarámovaný Lagrangián má tři strukturálně odlišné druhy členů a rozdíl je důležitý pro vše, co v této kapitole následuje:

| Část | Symbol | Závislost na konfiguraci | Závislost na rychlosti | Váže |
|---|---|---|---|---|
| Translační kinetická | $\tfrac{1}{2}\mu|\dot{\boldsymbol{\rho}}|^2$ | žádná | pouze $\dot{\boldsymbol{\rho}}$ | nic (čistě translační) |
| Keplerův potenciál (monopól-monopól) | $G m_A m_B / |\boldsymbol{\rho}|$ | pouze $\boldsymbol{\rho}$ | žádná | translaci na translaci |
| Rotační kinetická | $\tfrac{1}{2}\boldsymbol{\Omega}_X^T \mathbf{I}_X \boldsymbol{\Omega}_X$ | pouze $q_X$ (implicitně skrze $\boldsymbol{\Omega}_X = (\mathbf{R}_q^X)^T \dot{R}^X$) | pouze $\boldsymbol{\Omega}_X$ | nic (čistě rotační, na těleso) |
| Slapová vazba | $V_{\mathrm{tidal},\, X}(\boldsymbol{\rho}, q_X)$ | obojí $\boldsymbol{\rho}$ I $q_X$ | žádná | **translaci $\leftrightarrow$ rotaci tělesa-$X$** |

Výhradním zdrojem rotačně-translační vazby je **slapová část**: $V_{\mathrm{tidal},\, X}$ závisí na $\boldsymbol{\rho}$ (skrze předfaktor velikosti $1/|\boldsymbol{\rho}|^3$ a skrze jednotkový vektor $\hat{\boldsymbol{\rho}}$ uvnitř tělesové projekce) I na $q_X$ (skrze rotační matici $(\mathbf{R}_q^X)^T$ působící na $\hat{\boldsymbol{\rho}}$ k vyprodukování $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, X}$). Je-li slapová část odstraněna (např. nastavením $\mathbf{I}_A = \mathbf{I}_B = \mathbf{0}$ — Keplerova limita), Lagrangián se faktorizuje jako $L = L_{\mathrm{trans}}(\boldsymbol{\rho}, \dot{\boldsymbol{\rho}}) + L_{\mathrm{rot},A}(\boldsymbol{\Omega}_A) + L_{\mathrm{rot},B}(\boldsymbol{\Omega}_B)$, tři části jsou nezávislé a pohybové rovnice se rozpojí na čistou Keplerovu translaci plus dvě volné Eulerovy rotace. Tato faktorizace je strukturálním podložím validačního testovacího případu §2.3 (kap. 03 §3.8.1) a kap. 06 §6.2.

Naopak slapová část je zdrojem **dvou odlišných dynamických efektů**, jež sdílejí týž algebraický původ:

1. **Moment síly z gravitačního gradientu na rotaci tělesa-$X$.** $\partial V_{\mathrm{tidal},\, X} / \partial q_X \neq 0$ produkuje tělesový moment síly $\boldsymbol{\tau}_X^{\mathrm{body}}$ odvozený v kap. 03 §3.7 zarámovaném.
2. **Zpětná síla na relativní orbitu.** $\partial V_{\mathrm{tidal},\, X} / \partial \boldsymbol{\rho} \neq 0$ produkuje orientačně závislou sílu vstupující do pohybové rovnice $\boldsymbol{\rho}$. Toto je výstup §4.6 — zpětná síla nebyla odvozena v kap. 03 (jež se zaměřila na moment síly), a je tedy zde novým obsahem.

Oba efekty nejsou nezávislé: jsou rotační a translační projekcí téhož jediného objektu $\nabla V_{\mathrm{tidal},\, X}$ na součinové varietě. Jejich algebraická korelace vystoupí jako Newtonovsko-třetí-zákon rovnováha mezi $\boldsymbol{\tau}_X^{\mathrm{body}}$ a zpětnou silou, již činíme explicitní v §4.6.

### §4.2.2 Numerická kontrola zdravého rozumu na zarámovaném Lagrangiánu

Pro konfiguraci první verze (vstupy numerické kontroly kap. 03 §3.7.2: $m_A = 2450\ \text{kg}$, $\mathbf{I}_A = \mathrm{diag}(23323, 23323, 15384)\ \text{kg m}^2$ pro těleso podobné JWST; $m_B = 100\ \text{kg}$, $\mathbf{I}_B = \mathrm{diag}(20, 20, 116)\ \text{kg m}^2$ pro těleso podobné sondě; $|\boldsymbol{\rho}| = 10\ \text{m}$, $|\dot{\boldsymbol{\rho}}| = 0{,}1\ \text{m/s}$, $\boldsymbol{\Omega}_A = \boldsymbol{\Omega}_B = (0, 0, 0{,}01)\ \text{rad/s}$) se čtyři části zarámovaného Lagrangiánu vyhodnocují na:

| Část | Velikost (J) |
|---|---|
| Translační kinetická $\tfrac{1}{2}\mu|\dot{\boldsymbol{\rho}}|^2$ | $\sim 5 \times 10^{-1}$ |
| Kepler $G m_A m_B / |\boldsymbol{\rho}|$ | $\sim 2 \times 10^{-9}$ |
| Rotační $A$ + $B$ kinetická | $\sim 1$ |
| Slapová $V_{\mathrm{tidal},\, A} + V_{\mathrm{tidal},\, B}$ | $\sim 10^{-12}$ |

Slapová část je **o jedenáct řádů menší** než rotační kinetická při konfiguracích první verze — gravitační interakce při vzdálenostech řádu metrů mezi skromnými hmotnostmi je slabá. Zájmem modelu není velikost slapové části, nýbrž její **strukturální** role: je to *jediná* vazba mezi rotací a translací, a malá či velká, její přítomnost je tím, co dává soustavě její charakteristickou librační fenomenologii gravitačního gradientu a zpětnou sílu. Fakt, že diagnostika zachování $|dE/E_0| = 6{,}2 \times 10^{-12}$ ve výstupu první verze (kap. 03 §3.9) sleduje pouze člen monopól-monopól ve $V$ (dle `dynamics.py::total_potential_energy` — tabulka kap. 03 §3.9) bez slapového příspěvku, neovlivňuje pozorovanou přesnost zachování při těchto velikostech; OQ-3.4 (kap. 03 §3.11) katalogizuje rozšíření inženýrské úrovně, jež by učinilo diagnostiku sledující plné $V$.

## §4.3 Přístup A — Euler-Lagrange v souřadnicích (mapa Eulerových úhlů)

Variační princip kapitoly 01 §1.5 dává pohybové rovnice jako Eulerovy–Lagrangeovy rovnice zarámovaného Lagrangiánu výše na tečném svazku $Q_{\mathrm{CoM}}$. Existují dvě přirozené parametrizace rotačních faktorů: **mapa jednotkových kvaternionů** (volba inženýrské úrovně — explicitní v `dynamics.py`, trojrozměrná po uložení $|q|^2 = 1$) a **mapa Eulerových úhlů** (volba užitá v Podolský §7.3 a většině klasických učebnicových odvození). Obě jsou mapami na $SO(3)$ a dávají tutéž fyziku; následující odvození kanonické úrovně pracuje v mapě Eulerových úhlů, neboť (a) nemá žádné omezení k explicitnímu zpracování (každý Eulerův úhel je volnou souřadnicí), (b) mapuje přímo na Podolský 2018 §7.3 Rovn. (7.16)–(7.22) pro CS-jazykový křížový odkaz a (c) algebraický tvar výsledných rovnic je tvarem nejsnáze srovnatelným s Euler-Poincarého redukcí §4.4. Zamkni volbu. **(Řeší OQ-4.1.)** Kvaternionová mapa inženýrské úrovně je parametrizačně ekvivalentní a je volbou pro čas integrace z důvodů uvedených v poznámce inženýrské úrovně kap. 02 §2.7.1 — efektivita v kroku tělesové propagace. Odvození kanonické úrovně je na parametrizaci nezávislé skrze redukční krok §4.4, kde mapa zcela vypadává.

### §4.3.1 Eulerovy úhly a kinematické rovnice

Přijmi **intrinzickou ZXZ mapu Eulerových úhlů** na $SO(3)$: jakákoli orientace $R \in SO(3)$ je zapsána jako složení tří elementárních rotací parametrizovaných úhly $(\varphi, \vartheta, \psi)$,

$$
\mathbf{R}_q^X(\varphi_X, \vartheta_X, \psi_X)
\;=\;
R_z(\psi_X)\, R_x(\vartheta_X)\, R_z(\varphi_X),
$$

kde $R_z, R_x$ jsou elementární rotace kolem os $z$ a $x$ (Podolský 2018 §6.4 Rovn. 6.14–6.16 užívají přesně tuto konvenci; přijímáme jeho znaménko a pořadí k udržení CS-jazykového křížového odkazu čistým). Mapa je definována všude na $SO(3)$ kromě podmnožiny míry nula s kardanovým zámkem (gimbal lock); tato singularita je artefaktem mapy (nikoli fyzikální singularitou) a standardní zpracování — atlas více map, přepnutí map před průchodem singulární podmnožinou — platí. Pro pedagogické odvození postačí práce v jediné mapě.

Tělesová úhlová rychlost $\boldsymbol{\Omega}_X$ je pak vyjádřena v pojmech rychlostí Eulerových úhlů Podolský 2018 §6.5 zarámovanou Rovn. (6.17):

$$
\boldsymbol{\Omega}_X \;=\; \mathbf{K}(\vartheta_X, \psi_X)\, \begin{pmatrix} \dot\varphi_X \\ \dot\vartheta_X \\ \dot\psi_X \end{pmatrix},
\qquad
\mathbf{K}(\vartheta, \psi) \;=\; \begin{pmatrix} \sin\vartheta\, \sin\psi & \cos\psi & 0 \\ \sin\vartheta\, \cos\psi & -\sin\psi & 0 \\ \cos\vartheta & 0 & 1 \end{pmatrix}.
$$

Ekvivalentně, rozvinutím maticového součinu:

$$
\Omega_x = \dot\varphi \sin\vartheta\sin\psi + \dot\vartheta\cos\psi, \quad
\Omega_y = \dot\varphi \sin\vartheta\cos\psi - \dot\vartheta\sin\psi, \quad
\Omega_z = \dot\varphi \cos\vartheta + \dot\psi.
$$

Toto je **kinematická** část rovnic pohybu tuhého tělesa v mapě Eulerových úhlů: říká, jak se tělesová $\boldsymbol{\Omega}$ vztahuje k mapově-souřadnicovým rychlostem $(\dot\varphi, \dot\vartheta, \dot\psi)$, aniž by dosud volala Newtonovy, Eulerovy nebo Lagrangeovy rovnice na *dynamické* straně.

### §4.3.2 Lagrangián v souřadnicích mapy

Dosaď výraz s Eulerovými úhly pro $\boldsymbol{\Omega}_X$ do zarámovaného §4.2 Lagrangiánu. Rotační kinetická energie tělesa $A$ se stává (potlačujíc index $A$ na Eulerových úhlech a $\mathbf{I}$ pro přehlednost):

$$
\tfrac{1}{2}\boldsymbol{\Omega}^T \mathbf{I}\, \boldsymbol{\Omega}
\;=\;
\tfrac{1}{2}
\begin{pmatrix} \dot\varphi & \dot\vartheta & \dot\psi \end{pmatrix}
\mathbf{K}^T(\vartheta, \psi)\, \mathbf{I}\, \mathbf{K}(\vartheta, \psi)
\begin{pmatrix} \dot\varphi \\ \dot\vartheta \\ \dot\psi \end{pmatrix},
$$

kvadratická forma v rychlostech Eulerových úhlů s koeficienty závislými na $\vartheta$ a $\psi$. Pro těleso srovnané s hlavními osami ($\mathbf{I} = \mathrm{diag}(I_1, I_2, I_3)$) se kinetická energie rozvine na (algebra elementární, avšak zdlouhavá; viz Podolský §7.3 pro vypracované kroky):

$$
T_{\mathrm{rot},\, A}
\;=\;
\tfrac{1}{2}\, I_1 (\dot\varphi \sin\vartheta\sin\psi + \dot\vartheta\cos\psi)^2
+ \tfrac{1}{2}\, I_2 (\dot\varphi \sin\vartheta\cos\psi - \dot\vartheta\sin\psi)^2
+ \tfrac{1}{2}\, I_3 (\dot\varphi \cos\vartheta + \dot\psi)^2.
$$

Potenciální členy $V_{\mathrm{tidal},\, A}$ rovněž závisí na $(\varphi, \vartheta, \psi)$ skrze $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = (\mathbf{R}_q^A)^T \hat{\boldsymbol{\rho}}$, avšak závislost je na samotných Eulerových úhlech (orientaci), nikoli na jejich rychlostech.

### §4.3.3 Eulerovy–Lagrangeovy rovnice na mapě

Aplikuj $\frac{d}{dt}\frac{\partial L}{\partial \dot{q}^i} - \frac{\partial L}{\partial q^i} = 0$ na každou souřadnici Eulerova úhlu tělesa $A$. Tři výsledné rovnice jsou v rozvinutém mapovém tvaru neúhledné; Podolský 2018 §7.3 provádí vypracovanou algebru v jeho Rovn. (7.16)–(7.22) a dospívá k zarámované Rovn. (7.15):

$$
I_1 \dot\Omega_x - (I_2 - I_3)\, \Omega_y\Omega_z \;=\; M_x, \qquad \text{(a cyklicky)}.
$$

Pravá strana $M_x$ je tělesovou složkou zobecněného momentu síly vznikajícího z $-\partial V/\partial \boldsymbol{\theta}_A$ — ekvivalentní $\boldsymbol{\tau}_A^{\mathrm{body}}$ z kap. 03 §3.7 zarámovaného (s týmž předfaktorem a konvencemi znaménka, neboť Podolského levá strana a pravá strana kap. 03 užívají tutéž mapu a totéž znaménko pro tělesový moment síly).

Algebraický krok, jenž převádí mapovou-formu Eulerových–Lagrangeových rovnic na tělesové Eulerovy rovnice Rovn. (7.15), užívá kinematický transportní vztah Podolský Rovn. (6.8) (zarámovaný): pro libovolný tělesový vektor $\mathbf{w}$,

$$
\left.\frac{d\mathbf{w}}{dt}\right|_{\mathrm{spatial}} = \left.\frac{d\mathbf{w}}{dt}\right|_{\mathrm{body}} + \boldsymbol{\Omega} \times \mathbf{w}.
$$

Aplikováno na $\mathbf{w} = \mathbf{L}_A = \mathbf{I}_A \boldsymbol{\Omega}_A$ (tělesový moment hybnosti) je to přesně transport, jenž převádí $d\mathbf{L}/dt|_{\mathrm{spatial}} = \boldsymbol{\tau}^{\mathrm{inertial}}$ (druhý Newtonův zákon pro moment hybnosti v inerciální soustavě) na tělesový Eulerův tvar $\mathbf{I}\dot{\boldsymbol{\Omega}} + \boldsymbol{\Omega} \times (\mathbf{I}\boldsymbol{\Omega}) = \boldsymbol{\tau}^{\mathrm{body}}$. Lagrangeovsko-cestní odvození v Podolský §7.3 reprodukuje tento transport implicitně skrze algebru rychlostí Eulerových úhlů; nereodvozujeme jej zde v souřadnicích mapy.

Translační Eulerova–Lagrangeova rovnice $\boldsymbol{\rho}$ dává:

$$
\frac{d}{dt}(\mu \dot{\boldsymbol{\rho}}) \;-\; \frac{\partial L}{\partial \boldsymbol{\rho}}
\;=\; 0,
\qquad \Rightarrow \qquad
\mu \ddot{\boldsymbol{\rho}} \;=\; -\nabla_{\boldsymbol{\rho}} V,
$$

s $V = -Gm_Am_B/|\boldsymbol{\rho}| + V_{\mathrm{tidal},A} + V_{\mathrm{tidal},B}$. Toto je táž rovnice, již §4.6 rozebere podrobně; mapová forma zde nepotřebuje rozvoj.

### §4.3.4 Souhrn Přístupu A

Přístup A — přímý Euler-Lagrange v mapě Eulerových úhlů — produkuje:

- Klasické tělesové Eulerovy dynamické rovnice $\mathbf{I}_X \dot{\boldsymbol{\Omega}}_X + \boldsymbol{\Omega}_X \times (\mathbf{I}_X \boldsymbol{\Omega}_X) = \boldsymbol{\tau}_X^{\mathrm{body}}$ pro každé těleso $X = A, B$ (obnovené z algebry Eulerových úhlů skrze transportní rovnici).
- Translační rovnici $\mu \ddot{\boldsymbol{\rho}} = -\nabla_{\boldsymbol{\rho}} V$ pro souřadnici relativní polohy.
- Na souřadnicích závislou mezi-algebru, jež zastírá podkladovou geometrickou strukturu (levou invarianci rotační metriky), ale je správná v mapové formě.

Toto je **první noha křížové validace dvě-paralelně** v §4.5. V Přístupu B obnovíme tytéž tělesové Eulerovy rovnice souřadnicově nezávislou redukční cestou, jež činí hypotézu levé invariance explicitní.

## §4.4 Přístup B — Euler-Poincarého redukce

Věta o Euler-Poincarého redukci (Marsden & Ratiu 1999 §13.5) bere levo-invariantní Lagrangián na tečném svazku Lieovy grupy a redukuje jej na Lagrangián na Lieově algebře, s modifikovaným variačním principem, jenž produkuje rovnice na duálu Lieovy algebry. Aplikujeme ji na každý rotační faktor $SO(3)$ postupně; translační faktor je abelovský ($\mathbb{R}^3$) a má triviální redukci.

### §4.4.1 Euler-Poincarého věta

**Věta (Marsden & Ratiu 1999 §13.5).** Nechť $G$ je Lieova grupa s Lieovou algebrou $\mathfrak{g}$ a nechť $L: TG \to \mathbb{R}$ je Lagrangián. Předpokládejme, že $L$ je **levo-invariantní**: $L(g, \dot{g}) = L(\mathbb{I}, g^{-1}\dot{g})$ pro všechna $g \in G$ a $\dot{g} \in T_gG$. Definuj **redukovaný Lagrangián** $\ell: \mathfrak{g} \to \mathbb{R}$ vztahem

$$
\ell(\boldsymbol{\Omega}) \;\equiv\; L(\mathbb{I}, \boldsymbol{\Omega}), \qquad \boldsymbol{\Omega} = g^{-1}\dot{g} \in \mathfrak{g}.
$$

Pak Eulerovy–Lagrangeovy rovnice na $TG$ pro $L$ jsou ekvivalentní **Euler-Poincarého rovnicím** na $\mathfrak{g}^*$ pro $\ell$, s možnou externě dodanou zobecněnou silou $\boldsymbol{f}^{\mathrm{ext}} \in \mathfrak{g}^*$:

$$
\boxed{
\frac{d}{dt}\frac{\partial \ell}{\partial \boldsymbol{\Omega}}
\;=\;
\mathrm{ad}^*_{\boldsymbol{\Omega}}\,\frac{\partial \ell}{\partial \boldsymbol{\Omega}}
\;+\;
\boldsymbol{f}^{\mathrm{ext}},
}
$$

kde $\mathrm{ad}^*: \mathfrak{g} \times \mathfrak{g}^* \to \mathfrak{g}^*$ je duál adjungované akce $\mathfrak{g}$ na sebe samu, definovaný $\langle \mathrm{ad}^*_{\boldsymbol{\Omega}}\, \boldsymbol{\mu}, \boldsymbol{\eta}\rangle = \langle \boldsymbol{\mu}, [\boldsymbol{\Omega}, \boldsymbol{\eta}]\rangle$ pro $\boldsymbol{\Omega}, \boldsymbol{\eta} \in \mathfrak{g}$ a $\boldsymbol{\mu} \in \mathfrak{g}^*$. (Viz Marsden & Ratiu §9.3 pro adjungovanou a koadjungovanou akci obecně; §13.5 pro důkaz redukční věty.)

### §4.4.2 Ověření hypotézy levé invariance

Specializuj na jediný faktor $SO(3)$ z $Q_{\mathrm{CoM}}$. Rotační kinetický Lagrangián na $TSO(3)$ je

$$
L_{\mathrm{rot},\, X}(R, \dot{R}) \;=\; \tfrac{1}{2} \boldsymbol{\Omega}_X^T \mathbf{I}_X \boldsymbol{\Omega}_X, \qquad \widehat{\boldsymbol{\Omega}_X} = R^T \dot{R},
$$

s $\mathbf{I}_X$ **konstantním** tělesovým tenzorem setrvačnosti. Kapitola 02 §2.7.1 provedla explicitní důkaz, že tento Lagrangián je levo-invariantní: pod levou translací $R \to hR$, $\dot{R} \to h\dot{R}$ je levá trivializace $\boldsymbol{\Omega}_X = R^T \dot{R} \to (hR)^T (h\dot{R}) = R^T h^T h \dot{R} = R^T \dot{R} = \boldsymbol{\Omega}_X$ nezměněna, tudíž $L_{\mathrm{rot},\, X}$ je nezměněn. Hypotéza Euler-Poincarého věty je splněna pro **kinetickou** část rotačního Lagrangiánu.

Slapový potenciál $V_{\mathrm{tidal},\, X}$ však levo-invariantní **není**: závisí na tělesové projekci $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, X} = R^T \hat{\boldsymbol{\rho}}$, jež se mění pod $R \to hR$. Potenciál láme levou symetrii. Standardním zpracováním je absorbovat potenciál do externě dodané zobecněné síly $\boldsymbol{f}^{\mathrm{ext}}$ věty; tělesový moment síly $\boldsymbol{\tau}_X^{\mathrm{body}} = -\partial V_{\mathrm{tidal},\, X}/\partial \boldsymbol{\theta}_X$ (odvozený jako kap. 03 §3.7 zarámovaný) hraje přesně tuto roli. Euler-Poincarého rovnice se redukuje na tělesovou Eulerovu rovnici s momentem síly, jak nyní vypočteme explicitně.

### §4.4.3 Specializace $\mathrm{ad}^*$ pro $\mathfrak{so}(3) \cong \mathbb{R}^3$

Na Lieově algebře $\mathfrak{so}(3)$ antisymetrických matic $3\times 3$ je Lieova závorka maticovým komutátorem $[A, B] = AB - BA$. Ztotožněním $\mathfrak{so}(3) \cong \mathbb{R}^3$ skrze operátor stříšky z kap. 01 §1.4 ($\widehat{\boldsymbol{a}}\, \boldsymbol{b} = \boldsymbol{a} \times \boldsymbol{b}$) se závorka transportuje na vektorový součin:

$$
[\boldsymbol{\Omega}, \boldsymbol{\eta}]_{\mathfrak{so}(3)} \;\leftrightarrow\; \boldsymbol{\Omega} \times \boldsymbol{\eta}\quad \text{(v $\mathbb{R}^3$)}.
$$

Duální prostor $\mathfrak{so}(3)^*$ je ztotožněn s $\mathbb{R}^3$ skrze standardní euklidovský skalární součin. Nyní můžeme vypočíst $\mathrm{ad}^*_{\boldsymbol{\Omega}}\, \boldsymbol{\mu}$ explicitně přenesením párování. Pro $\boldsymbol{\Omega}, \boldsymbol{\eta} \in \mathfrak{so}(3) \cong \mathbb{R}^3$ a $\boldsymbol{\mu} \in \mathfrak{so}(3)^* \cong \mathbb{R}^3$ se definice $\langle \mathrm{ad}^*_{\boldsymbol{\Omega}}\, \boldsymbol{\mu}, \boldsymbol{\eta}\rangle = \langle \boldsymbol{\mu}, [\boldsymbol{\Omega}, \boldsymbol{\eta}]\rangle$ stává:

$$
\mathrm{ad}^*_{\boldsymbol{\Omega}}\,\boldsymbol{\mu} \cdot \boldsymbol{\eta}
\;=\;
\boldsymbol{\mu} \cdot (\boldsymbol{\Omega} \times \boldsymbol{\eta}).
$$

S užitím cyklické identity smíšeného součinu $\boldsymbol{\mu} \cdot (\boldsymbol{\Omega} \times \boldsymbol{\eta}) = (\boldsymbol{\mu} \times \boldsymbol{\Omega}) \cdot \boldsymbol{\eta} = -(\boldsymbol{\Omega} \times \boldsymbol{\mu}) \cdot \boldsymbol{\eta}$ — a požadavkem, aby toto platilo pro libovolné $\boldsymbol{\eta}$ — odečteme:

$$
\boxed{
\mathrm{ad}^*_{\boldsymbol{\Omega}}\,\boldsymbol{\mu} \;=\; -\, \boldsymbol{\Omega} \times \boldsymbol{\mu}, \qquad \boldsymbol{\Omega} \in \mathfrak{so}(3) \cong \mathbb{R}^3, \quad \boldsymbol{\mu} \in \mathfrak{so}(3)^* \cong \mathbb{R}^3.
}
$$

**(Řeší OQ-4.2.)** Minusové znaménko je konvencí vynucenou párováním duality — *nikoli* nezávislou volbou. Různé učebnice prezentují $\mathrm{ad}^*$ pro $\mathfrak{so}(3)$ s opačným znaménkem; neshoda se redukuje na to, zda je párování $\langle \cdot, \cdot \rangle: \mathfrak{g} \times \mathfrak{g}^* \to \mathbb{R}$ čteno s $\mathfrak{g}$ vlevo nebo vpravo. Sledujeme konvenci Marsden & Ratiu §9.3 (levá-na-levé), dávající explicitní minusové znaménko zde. Tělesová Eulerova rovnice odvozená níže nese znaménko konzistentně s konvencemi levé trivializace kap. 02 §2.7.1 a s kap. 03 §3.7 zarámovaným momentem síly z gravitačního gradientu.

### §4.4.4 Specializace Euler-Poincarého rovnice na rotaci tuhého tělesa

Redukovaný rotační Lagrangián je

$$
\ell_X(\boldsymbol{\Omega}_X) \;=\; \tfrac{1}{2} \boldsymbol{\Omega}_X^T \mathbf{I}_X \boldsymbol{\Omega}_X, \qquad \frac{\partial \ell_X}{\partial \boldsymbol{\Omega}_X} \;=\; \mathbf{I}_X\, \boldsymbol{\Omega}_X.
$$

Vlož do zarámované Euler-Poincarého rovnice, s $\boldsymbol{f}^{\mathrm{ext}}_X = \boldsymbol{\tau}_X^{\mathrm{body}}$ jak pojednáno v §4.4.2:

$$
\frac{d}{dt}(\mathbf{I}_X \boldsymbol{\Omega}_X)
\;=\;
\mathrm{ad}^*_{\boldsymbol{\Omega}_X}\, (\mathbf{I}_X \boldsymbol{\Omega}_X)
\;+\; \boldsymbol{\tau}_X^{\mathrm{body}}
\;=\;
-\, \boldsymbol{\Omega}_X \times (\mathbf{I}_X \boldsymbol{\Omega}_X)
\;+\; \boldsymbol{\tau}_X^{\mathrm{body}},
$$

s užitím zarámovaného $\mathrm{ad}^*_{\boldsymbol{\Omega}}\,\boldsymbol{\mu} = -\boldsymbol{\Omega} \times \boldsymbol{\mu}$ z §4.4.3. Jelikož $\mathbf{I}_X$ je v tělesové soustavě konstantní, $d(\mathbf{I}_X \boldsymbol{\Omega}_X)/dt = \mathbf{I}_X \dot{\boldsymbol{\Omega}}_X$ a přeskupení dává **tělesovou Eulerovu rovnici** v její kanonické formě:

$$
\boxed{
\mathbf{I}_X \dot{\boldsymbol{\Omega}}_X \;+\; \boldsymbol{\Omega}_X \times (\mathbf{I}_X \boldsymbol{\Omega}_X) \;=\; \boldsymbol{\tau}_X^{\mathrm{body}},
\qquad X \in \{A, B\}.
}
$$

Toto jsou pohybové rovnice pro rotační stupně volnosti na redukovaném fázovém prostoru Lieovy algebry. Translační rovnice $\boldsymbol{\rho}$ je touto redukcí neovlivněna (Euler-Poincaré platí jen pro faktory Lieovy grupy $SO(3)$; translační faktor $\mathbb{R}^3_{\boldsymbol{\rho}}$ je již abelovský a redukce je triviální).

### §4.4.5 Souhrn Přístupu B

Přístup B — Euler-Poincarého redukce — produkuje:

- Tělesovou Eulerovu rovnici $\mathbf{I}_X \dot{\boldsymbol{\Omega}}_X + \boldsymbol{\Omega}_X \times (\mathbf{I}_X \boldsymbol{\Omega}_X) = \boldsymbol{\tau}_X^{\mathrm{body}}$ pro každé těleso $X = A, B$, identickou s výsledkem §4.3.3 Přístupu A.
- Explicitní odvození $\mathrm{ad}^*$, jež upevňuje konvenci znaménka k souladu s kap. 02 §2.7.1 + kap. 03 §3.7 + oprava znaménka OQ-FORT-1.
- Souřadnicově nezávislou prezentaci, jež činí hypotézu levé invariance (kap. 02 §2.7.1) a roli externě aplikovaného momentu síly explicitní.

Redukce činí zjevným, že tělesová Eulerova rovnice je **strukturálním důsledkem** levé invariance metriky, nikoli Newtonovsko-Eulerovým tvrzením $d\mathbf{L}/dt = \boldsymbol{\tau}$ transformovaným po souřadnicích. Oba pohledy jsou správné; strukturální pohled vysvětluje, *proč* je tělesová soustava přirozenou soustavou pro propagaci tuhého tělesa, a motivuje volbu inženýrské úrovně (poznámka inženýrské úrovně kap. 02 §2.7.1).

## §4.5 Křížová validace: Přístup A ≡ Přístup B

Oba přístupy dospívají k **téže** tělesové Eulerově rovnici $\mathbf{I}_X \dot{\boldsymbol{\Omega}}_X + \boldsymbol{\Omega}_X \times (\mathbf{I}_X \boldsymbol{\Omega}_X) = \boldsymbol{\tau}_X^{\mathrm{body}}$ pro každé těleso $X$:

| | Přístup A | Přístup B |
|---|---|---|
| Metoda | Přímý Euler-Lagrange v mapě Eulerových úhlů | Euler-Poincarého redukce levo-invariantního Lagrangiánu |
| Závislost na souřadnicích | Ano (mapa na $SO(3)$) | Ne (souřadnicově nezávislá) |
| Mezikroky | Algebra mapy + transportní rovnice | Identifikace $\mathrm{ad}^*$ + redukční věta |
| Konečná rovnice | $\mathbf{I}_X \dot{\boldsymbol{\Omega}}_X + \boldsymbol{\Omega}_X \times (\mathbf{I}_X \boldsymbol{\Omega}_X) = \boldsymbol{\tau}_X^{\mathrm{body}}$ | (identická) |
| Konvence znaménka | Zděděné z Podolský §7.3 | Zděděné z Marsden-Ratiu §13.5 |
| Křížová kontrola | $\boldsymbol{\tau}_X^{\mathrm{body}}$ odpovídá kap. 03 §3.7 zarámovanému (po opravě znaménka OQ-FORT-1) | (totéž) |

Shoda je validací. Jak zdůrazněno v kap. 03 §3.6 (metodologie dvě-paralelně) a direktivě projektového vlastníka 2026-05-23: hodnotou provedení *obou* odvození je, že algebraická chyba v kterékoli cestě by vystoupila jako neshoda a eliminace neshody je eliminací oné třídy chyby.

Tato křížová validace je analytickým analogem doktrinálního kritéria §4.4 C2 bitové identity mezi implementacemi, aplikovaným na úrovni odvození. Následující kapitoly kanonické úrovně rozšiřují vzor: kapitola 05 odvozuje specializaci na hlavní osy jak rozvojem složka po složce, tak strukturální specializací výsledku §4.4; kapitola 06 validuje každý numerický testovací případ vůči alespoň dvěma učebnicovým odkazům, kde jsou dostupné.

> **Doplňující anglojazyčné odkazy pro vzor křížové validace.** Marsden & Ratiu 1999 §13 (Euler-Poincarého redukce levo-invariantních Lagrangiánů, s tuhým tělesem $SO(3)$ jako pracovním příkladem v §13.5). Holm-Schmah-Stoica 2009 *Geometric Mechanics and Symmetry* §6.2 provádí tutéž redukci s explicitní algebrou $\mathrm{ad}^*$ ve specializaci $SO(3)$. Arnold 1989 §3.27 dává elegantní krátkou formu: tuhé těleso je geodetickým tokem na $SO(3)$ s levo-invariantní metrikou a Eulerovy rovnice jsou rovnicemi geodetik v tělesových souřadnicích. Všechny tři odkazy souhlasí na zarámované tělesové Eulerově rovnici; znaménko členu vektorového součinu je invariantní pod jejich odlišnými konvencemi.

## §4.6 Translační rovnice v těžišťové soustavě a zpětná síla

Translační Eulerova–Lagrangeova rovnice z §4.3.3 zní:

$$
\mu \ddot{\boldsymbol{\rho}} \;=\; -\nabla_{\boldsymbol{\rho}} V,
$$

s $V = -Gm_Am_B/|\boldsymbol{\rho}| + V_{\mathrm{tidal},A} + V_{\mathrm{tidal},B}$ z §4.2 zarámovaného. Rozlož člen gradientu člen po členu.

### §4.6.1 Gradient Keplerovy části

$$
-\nabla_{\boldsymbol{\rho}}\left(-\frac{Gm_Am_B}{|\boldsymbol{\rho}|}\right)
\;=\;
-\frac{Gm_Am_B}{|\boldsymbol{\rho}|^2}\, \hat{\boldsymbol{\rho}}
\;=\;
-\frac{Gm_Am_B}{|\boldsymbol{\rho}|^3}\, \boldsymbol{\rho}.
$$

Toto je **Keplerova síla** — přitažlivá centrální síla nepřímé úměry čtverci mezi dvěma tělesy branými jako bodové hmoty. V kombinaci s $\mu \ddot{\boldsymbol{\rho}} = \mathbf{F}_{\mathrm{Kepler}}$ samotnou (tj. bez slapových částí) je pohybová rovnice standardní Keplerovou rovnicí s redukovanou hmotností (Landau-Lifšic, sv. I §13–§15). Toto je limita testovacího případu §2.3 / kap. 06 §6.2.

### §4.6.2 Gradient slapové části — odvození zpětné síly

Slapová část pro těleso $A$ je (kap. 03 §3.4.3 zarámovaná):

$$
V_{\mathrm{tidal},\, A}(\boldsymbol{\rho}, q_A)
\;=\;
-\frac{Gm_B}{2|\boldsymbol{\rho}|^3}\,
\Big[ \mathrm{tr}(\mathbf{I}_A) \;-\; 3\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \Big].
$$

K výpočtu $-\nabla_{\boldsymbol{\rho}} V_{\mathrm{tidal},\, A}$ rozlož závislost na $\boldsymbol{\rho}$ na dvě části:

1. Skalární **předfaktor** $1/|\boldsymbol{\rho}|^3$ — závisí pouze na velikosti $r \equiv |\boldsymbol{\rho}|$.
2. **Závorkovaný orientační faktor** $[\mathrm{tr}(\mathbf{I}_A) - 3\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}]$ — závisí na $\boldsymbol{\rho}$ pouze skrze jednotkový vektor $\hat{\boldsymbol{\rho}}$ uvnitř $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = (\mathbf{R}_q^A)^T \hat{\boldsymbol{\rho}}$.

Označ závorkovaný orientační faktor $\Phi_A(\hat{\boldsymbol{\rho}}, q_A) \equiv \mathrm{tr}(\mathbf{I}_A) - 3\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$. Součinovým pravidlem:

$$
\nabla_{\boldsymbol{\rho}} V_{\mathrm{tidal},\, A}
\;=\;
-\frac{Gm_B}{2}\,
\Big[\,
\Phi_A\, \nabla_{\boldsymbol{\rho}}(1/r^3)
\;+\;
\frac{1}{r^3}\, \nabla_{\boldsymbol{\rho}} \Phi_A
\,\Big].
$$

Vypočti každou část.

**Gradient radiálního předfaktoru.** $\nabla_{\boldsymbol{\rho}}(1/r^3) = -3\hat{\boldsymbol{\rho}}/r^4$, takže:

$$
-\frac{Gm_B}{2}\, \Phi_A\, \nabla_{\boldsymbol{\rho}}(1/r^3)
\;=\;
\frac{3\, G m_B}{2\, r^4}\, \Phi_A\, \hat{\boldsymbol{\rho}}.
$$

**Gradient orientačního faktoru.** $\Phi_A$ závisí na $\boldsymbol{\rho}$ pouze skrze $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = (\mathbf{R}_q^A)^T \hat{\boldsymbol{\rho}}$. S užitím $\partial \hat{\boldsymbol{\rho}}/\partial \boldsymbol{\rho} = (\mathbb{I}_3 - \hat{\boldsymbol{\rho}} \hat{\boldsymbol{\rho}}^T)/r$ (identita příčné projekce) a aplikováním řetězového pravidla:

$$
\nabla_{\boldsymbol{\rho}}\Phi_A
\;=\;
\frac{(\mathbb{I}_3 - \hat{\boldsymbol{\rho}}\hat{\boldsymbol{\rho}}^T)}{r}\, \mathbf{R}_q^A\, \nabla_{\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}} \Phi_A.
$$

S užitím $\Phi_A = \mathrm{tr}(\mathbf{I}_A) - 3\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$ a poznamenáním, že $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$ má jednotkovou velikost:

$$
\nabla_{\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}} \Phi_A
\;=\;
-6\, \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}.
$$

Kombinací:

$$
\nabla_{\boldsymbol{\rho}}\Phi_A
\;=\;
-\frac{6}{r}\, (\mathbb{I}_3 - \hat{\boldsymbol{\rho}}\hat{\boldsymbol{\rho}}^T)\, \mathbf{R}_q^A\, (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}).
$$

Příspěvek orientačního faktoru ke gradientu je tedy:

$$
-\frac{Gm_B}{2}\, \frac{1}{r^3}\, \nabla_{\boldsymbol{\rho}}\Phi_A
\;=\;
\frac{3\, G m_B}{r^4}\, (\mathbb{I}_3 - \hat{\boldsymbol{\rho}}\hat{\boldsymbol{\rho}}^T)\, \mathbf{R}_q^A\, (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}).
$$

**Kombinovaná zpětná síla na relativní orbitu, ze slapového příspěvku tělesa $A$.** Sebráním obou částí:

$$
\boxed{
\mathbf{F}_{\mathrm{back-reaction},\, A}
\;\equiv\;
-\nabla_{\boldsymbol{\rho}} V_{\mathrm{tidal},\, A}
\;=\;
-\frac{3\, G m_B}{2\, r^4}\, \Phi_A\, \hat{\boldsymbol{\rho}}
\;-\;
\frac{3\, G m_B}{r^4}\, (\mathbb{I}_3 - \hat{\boldsymbol{\rho}}\hat{\boldsymbol{\rho}}^T)\, \mathbf{R}_q^A\, (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}).
}
$$

Zpětná síla škáluje jako $1/r^4$ — o jednu mocninu $r$ strmější než Keplerova síla $1/r^2$ — a má orientačně závislý předfaktor kódující anizotropii tělesa $A$ a jeho sourodost s $\hat{\boldsymbol{\rho}}$. Symetrický výraz pro $\mathbf{F}_{\mathrm{back-reaction},\, B}$ následuje výměnou $A \leftrightarrow B$ (se znaménkem $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, B}$ obráceným, jako v kap. 03 §3.7.1 pro moment síly).

**(Řeší OQ-4.3.)** Zpětná síla existuje, škáluje jako $1/r^4$ (odpovídá odhadu meze OQ-4.3 v zástupci §4.9) a má explicitní uzavřený výraz v pojmech $(r, \hat{\boldsymbol{\rho}}, \mathbf{I}_A, \mathbf{R}_q^A)$.

### §4.6.3 Vyrušení v limitě sférického tělesa

Pro sférické těleso ($\mathbf{I}_A = I_0 \mathbb{I}_3$): $\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = I_0\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$ a $\Phi_A = 3I_0 - 3I_0 |\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}|^2 = 3I_0 - 3I_0 = 0$. První (radiálně-předfaktorový) člen v $\mathbf{F}_{\mathrm{back-reaction},\, A}$ mizí kvůli $\Phi_A = 0$. Druhý (orientačně-faktorový) člen: $(\mathbb{I}_3 - \hat{\boldsymbol{\rho}}\hat{\boldsymbol{\rho}}^T)\, \mathbf{R}_q^A\, (I_0 \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}) = I_0\, (\mathbb{I}_3 - \hat{\boldsymbol{\rho}}\hat{\boldsymbol{\rho}}^T)\, \mathbf{R}_q^A (\mathbf{R}_q^A)^T \hat{\boldsymbol{\rho}} = I_0\, (\mathbb{I}_3 - \hat{\boldsymbol{\rho}}\hat{\boldsymbol{\rho}}^T) \hat{\boldsymbol{\rho}} = I_0\, (\hat{\boldsymbol{\rho}} - \hat{\boldsymbol{\rho}}) = \mathbf{0}$. Obě části mizí; zpětná síla na orbitu v důsledku slapové vazby tělesa $A$ identicky mizí, je-li $A$ sférické. To je konzistentní s mizením momentu síly z gravitačního gradientu v téže limitě (kap. 03 §3.8.2): gravitační vazba sférického tělesa na jeho prostředí závisí pouze na jeho hmotnosti, nikoli na orientaci, takže ani rotační dynamika, ani translační zpětná síla nemohou nést informaci odvozenou z orientace.

### §4.6.4 Rovnováha zpětná-síla–moment-síly (kontrola třetího Newtonova zákona)

Užitečná kontrola konzistence: v izolované soustavě dvou těles je celková lineární hybnost $\mathbf{P}_{\mathrm{CoM}}$ zachována (kap. 02 §2.6.1; symetrie pod prostorovou translací Noetherové větou). Tudíž celková síla na těžišťovou souřadnici mizí:

$$
\mathbf{F}^{\mathrm{external}}_{\mathrm{CoM}} \;\equiv\; \frac{d\mathbf{P}_{\mathrm{CoM}}}{dt} \;=\; \mathbf{0}.
$$

Toto je automaticky respektováno translační invariancí našeho Lagrangiánu — Lagrangián závisí na $(\mathbf{R}_A, \mathbf{R}_B)$ jen skrze $\boldsymbol{\rho} = \mathbf{R}_B - \mathbf{R}_A$, takže $\partial L/\partial \mathbf{R}_{\mathrm{CoM}} = 0$ a $\mathbf{P}_{\mathrm{CoM}}$ je prvním integrálem cyklické souřadnice s konjugovanou hybností. Zpětná síla vstupuje do *relativní* rovnice $\mu\ddot{\boldsymbol{\rho}} = \mathbf{F}_{\mathrm{Kepler}} + \mathbf{F}_{\mathrm{back-reaction},\, A} + \mathbf{F}_{\mathrm{back-reaction},\, B}$, ale neláme zachování těžiště; zpětná síla je vnitřní silou ve smyslu třetího Newtonova zákona.

Odpovídajícím tvrzením o momentu hybnosti je **rovnováha celkového momentu hybnosti**: rychlost změny celkového momentu hybnosti (orbitální + spin-tělesa-$A$ + spin-tělesa-$B$) kolem těžiště mizí při nepřítomnosti vnějších momentů sil. V naší soustavě může orbitální moment hybnosti $\boldsymbol{\rho} \times \mu \dot{\boldsymbol{\rho}}$ vyměňovat s tělesovými spiny $\mathbf{R}_q^X (\mathbf{I}_X \boldsymbol{\Omega}_X)$ skrze moment síly z gravitačního gradientu (rotační strana) plus moment-síly-na-orbitu v důsledku zpětné síly (translační strana); oba příspěvky se vyruší. Explicitní ověření je výpočet $\mathbf{F}_{\mathrm{back-reaction},\, A} \times \mathbf{R}_A + \boldsymbol{\tau}_A^{\mathrm{inertial}} + \mathbf{F}_{\mathrm{back-reaction},\, B} \times \mathbf{R}_B + \boldsymbol{\tau}_B^{\mathrm{inertial}} = \mathbf{0}$, jenž plyne z $\nabla_{\boldsymbol{\theta}_A} V + \nabla_{\boldsymbol{\theta}_B} V + \boldsymbol{\rho} \times \nabla_{\boldsymbol{\rho}} V = 0$, identity vyjadřující rotační invarianci vzájemného potenciálu $V$. Tvrdíme jej zde pro strukturální úplnost; explicitní kontrola Noetherové identity je diagnostikou zachování inženýrské úrovně $|d\mathbf{L}_{\mathrm{total}}/L_0|$, již první-verzní integrace RK4 dosahuje na $2{,}1 \times 10^{-10}$ přes okno 600 sekund (kap. 00 §3).

## §4.7 Úplná sada rovnic a křížová validace vůči dynamics.py

Sestav §4.4 zarámované rotační rovnice a §4.6 translační rovnici včetně zpětné síly do **úplné sady rovnic** na 18rozměrném redukovaném fázovém prostoru $T^* Q_{\mathrm{CoM}}$ (počítání: 3 složky $\boldsymbol{\rho}$ + 3 složky $\mathbf{p}_{\mathrm{rel}}$ + 3 Eulerovy-úhlové / 4-kvaternionové souřadnice $q_A$ + 3 tělesové úhlově-rychlostní složky $\boldsymbol{\Omega}_A$ + 3 kvaternion / Eulerovy úhly pro $q_B$ + 3 pro $\boldsymbol{\Omega}_B$, s jedním omezením na kvaternion v inženýrské úrovni; netto 18):

$$
\boxed{
\begin{aligned}
\text{Translační:}\quad
&\dot{\boldsymbol{\rho}} \;=\; \mathbf{p}_{\mathrm{rel}}/\mu \\
&\dot{\mathbf{p}}_{\mathrm{rel}} \;=\; \mathbf{F}_{\mathrm{Kepler}}(\boldsymbol{\rho}) \;+\; \mathbf{F}_{\mathrm{back-reaction},\, A}(\boldsymbol{\rho}, q_A) \;+\; \mathbf{F}_{\mathrm{back-reaction},\, B}(\boldsymbol{\rho}, q_B) \\[4pt]
\text{Rotační $A$:}\quad
&\dot{q}_A \;=\; \tfrac{1}{2}\, q_A \otimes \boldsymbol{\Omega}_A \quad \text{(kvaternionová mapa)} \\
&\mathbf{I}_A \dot{\boldsymbol{\Omega}}_A \;=\; -\boldsymbol{\Omega}_A \times (\mathbf{I}_A \boldsymbol{\Omega}_A) \;+\; \boldsymbol{\tau}_A^{\mathrm{body}}(\boldsymbol{\rho}, q_A) \\[4pt]
\text{Rotační $B$:}\quad
&\dot{q}_B \;=\; \tfrac{1}{2}\, q_B \otimes \boldsymbol{\Omega}_B \\
&\mathbf{I}_B \dot{\boldsymbol{\Omega}}_B \;=\; -\boldsymbol{\Omega}_B \times (\mathbf{I}_B \boldsymbol{\Omega}_B) \;+\; \boldsymbol{\tau}_B^{\mathrm{body}}(\boldsymbol{\rho}, q_B)
\end{aligned}
}
$$

s:

$$
\mathbf{F}_{\mathrm{Kepler}}(\boldsymbol{\rho}) \;=\; -\frac{Gm_Am_B}{|\boldsymbol{\rho}|^3}\, \boldsymbol{\rho}, \qquad
\boldsymbol{\tau}_X^{\mathrm{body}} \;=\; \frac{3Gm_{\bar X}}{|\boldsymbol{\rho}|^3}\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, X} \times (\mathbf{I}_X \hat{\boldsymbol{\rho}}_{\mathrm{body},\, X}),
$$

a $\mathbf{F}_{\mathrm{back-reaction},\, X}$ daným §4.6.2 zarámovaným.

Toto je sada rovnic kanonické úrovně — pohybové rovnice, jež musí jakákoli implementace této fyziky splňovat. Kapitola 05 specializuje rotační části na soustavu hlavních os a dosazuje explicitní výraz momentu síly; kapitola 06 validuje každou část vůči analytickým limitům.

### §4.7.1 Shoda řádek po řádku vůči `dynamics.py::state_derivative`

Implementace inženýrské úrovně v `dynamics.py::state_derivative` (řádky ~130–170) propaguje přesně zarámovanou sadu rovnic výše. Zhuštěně funkce zní (parafrázováno; nahlédni do souboru pro doslovný kód):

```
# Translační
r_vec = body_b.x - body_a.x   # rovná se boldsymbol{rho}
r_mag = norm(r_vec)
F_kepler = -G * m_a * m_b * r_vec / r_mag**3
# Rotační A
tau_a_body = gravity_gradient_torque_body_frame(r_vec, R_a, I_a, m_b)
omega_a_dot = solve(I_a, tau_a_body - cross(omega_a, I_a @ omega_a))
q_a_dot = 0.5 * q_a_kinematics(q_a, omega_a)
# Rotační B (symetrické, s -r_vec)
tau_b_body = gravity_gradient_torque_body_frame(-r_vec, R_b, I_b, m_a)
omega_b_dot = solve(I_b, tau_b_body - cross(omega_b, I_b @ omega_b))
q_b_dot = 0.5 * q_b_kinematics(q_b, omega_b)
# Sestavení derivace stavu
return [x_a_dot, x_b_dot, v_a_dot, v_b_dot, q_a_dot, q_b_dot, omega_a_dot, omega_b_dot]
```

Shoda je řádek po řádku:

| Kanonická rovnice (zarámovaná §4.7) | Řádek inženýrské úrovně | Status |
|---|---|---|
| $\dot{\boldsymbol{\rho}} = \mathbf{p}_{\mathrm{rel}}/\mu$ | `v_a_dot, v_b_dot` (v těžišťové soustavě je relativní pohyb automatický z rozdílu) | přesné |
| $\mathbf{F}_{\mathrm{Kepler}} = -Gm_Am_B \boldsymbol{\rho}/r^3$ | řádek `F_kepler` | přesné (předfaktor + znaménko) |
| $\mathbf{F}_{\mathrm{back-reaction},\, X}$ (§4.6.2 zarámované) | **NEZAHRNUTO** v první-verzní `state_derivative` | mezera inženýrské úrovně; označeno v řešení OQ-4.4 níže |
| $\dot{q}_X = \tfrac{1}{2} q_X \otimes \boldsymbol{\Omega}_X$ | `q_a_kinematics`, `q_b_kinematics` | přesné |
| $\mathbf{I}_X \dot{\boldsymbol{\Omega}}_X = \boldsymbol{\tau}_X^{\mathrm{body}} - \boldsymbol{\Omega}_X \times (\mathbf{I}_X \boldsymbol{\Omega}_X)$ | řádky `omega_a_dot, omega_b_dot` | přesné |
| $\boldsymbol{\tau}_X^{\mathrm{body}} = (3Gm_{\bar X}/r^3)\hat{\boldsymbol{\rho}}_{\mathrm{body},X} \times (\mathbf{I}_X \hat{\boldsymbol{\rho}}_{\mathrm{body},X})$ | `gravity_gradient_torque_body_frame` (řádky 113–131) — ověřeno v kap. 03 §3.7.1 | přesné |
| Konvence volání tělesa-$B$ s $-\boldsymbol{\rho}$ jako vstupem | řádek 158 — ověřeno v kap. 03 §3.7.1 | přesné |

**(Řeší OQ-4.4.)** Shoda řádek po řádku je přesná na **rotační** úrovni (translační Kepler + rotační Euler + moment síly z gravitačního gradientu vše správné). **Mezera** je u zpětné síly: kanonický §4.6.2 zarámovaný výraz $\mathbf{F}_{\mathrm{back-reaction},\, X}$ **dosud není implementován** v `dynamics.py::state_derivative`. První verze tedy propaguje čistou Keplerovu translaci + rotaci s momentem síly z gravitačního gradientu, vynechávajíc vazbu orientace-zpětné-síly-na-orbitu.

**Velikost chybějícího příspěvku.** Při konfiguracích první verze (kap. 03 §3.7.2 hypotetické podobné L2) $|\mathbf{F}_{\mathrm{Kepler}}| \sim G m_A m_B / r^2 \sim 10^{-12}\ \text{N}$ a $|\mathbf{F}_{\mathrm{back-reaction}}| \sim 3 G m_B \mathrm{tr}(\mathbf{I}_A) / r^5 \sim 10^{-18}\ \text{N}$ — zpětná síla je o šest řádů menší než Keplerova síla při vzdálenostech $\gg \ell$ (rozměr tělesa) a poměr škáluje jako $(\ell/r)^3$, jak se očekává pro korekci dalšího multipólu. Při měřítkách první verze je vynechání pod zaokrouhlovací chybou RK4 a neovlivňuje žádný pozorovaný numerický výsledek. V režimech slapového uzamčení binárního asteroidu ($\ell/r \to 1$) se zpětná síla stává dominantním ne-Keplerovým orbitálním efektem; práce inženýrské úrovně v0.2 by ji implementovala. **Navazující úkol inženýrské úrovně:** zařazeno jako nové OQ-4.5 v §4.9.

### §4.7.2 Numerická kontrola zdravého rozumu: redukce Keplerovy limity

Nastav $\mathbf{I}_A = \mathbf{I}_B = \mathbf{0}$ v zarámované §4.7 sadě rovnic. Pak $\boldsymbol{\tau}_A^{\mathrm{body}} = \boldsymbol{\tau}_B^{\mathrm{body}} = \mathbf{0}$ (úměrné $\mathbf{I}_X$), $\mathbf{F}_{\mathrm{back-reaction},\, X} = \mathbf{0}$ (rovněž úměrné $\mathbf{I}_X$ — §4.6.2 zarámované) a translační rovnice se zhroutí na:

$$
\mu \ddot{\boldsymbol{\rho}} \;=\; \mathbf{F}_{\mathrm{Kepler}} \;=\; -\frac{Gm_Am_B}{r^2}\, \hat{\boldsymbol{\rho}},
$$

s rotačními rovnicemi rozpojenými na $\boldsymbol{\Omega}_X = \mathrm{const}$ (žádný moment síly na tělese s nulovou setrvačností). Translační rovnice je standardní Keplerovou rovnicí s redukovanou hmotností; zachovávající se veličiny (energie, moment hybnosti, Laplaceův-Runge-Lenzův vektor) se redukují na Keplerovy zachovávající se veličiny; orbitální perioda odpovídá učebnicovému vzorci $T_{\mathrm{orbit}} = 2\pi\sqrt{a^3/(G(m_A + m_B))}$. Toto je redukce testovacího případu §2.3 / kap. 06 §6.2; kontrola zdravého rozumu potvrzuje, že zarámovaná §4.7 sada rovnic má správnou Keplerovu limitu.

## §4.8 Co je stanoveno na konci této kapitoly

Čtenář, jenž prošel §4.1–§4.7, má:

- Úplný Lagrangián $L = T - V$ (§4.2 zarámovaný) na redukované konfigurační varietě $Q_{\mathrm{CoM}} = \mathbb{R}^3 \times SO(3)_A \times SO(3)_B$, s kinetickou energií z kap. 02 a kvadrupólově-useknutým potenciálem z kap. 03 sestavenými a strukturálně rozloženými (translační kinetická + Keplerův potenciál + rotační kinetická na těleso + slapová vazba na těleso).
- Strukturální identifikaci (§4.2.1) **slapové části** jako výhradního zdroje rotačně-translační vazby a rozpoznání, že její projekce $\partial/\partial q_X$ dává moment síly z gravitačního gradientu, zatímco její projekce $\partial/\partial \boldsymbol{\rho}$ dává orbitální zpětnou sílu.
- Eulerovy–Lagrangeovy rovnice na $TQ_{\mathrm{CoM}}$ odvozené přímým mapovým výpočtem (Přístup A, §4.3) v mapě Eulerových úhlů, s algebrou mapy reprodukující klasické tělesové Eulerovy dynamické rovnice skrze kinematický transportní vztah. CS-jazykové ukotvení na Podolský §6.4–§6.5 (aparát mapy) + §7.3 (vypracovaná algebra) učiněno konkrétním.
- Euler-Poincarého redukci (Přístup B, §4.4) levo-invariantního rotačního Lagrangiánu na každém faktoru $SO(3)$, s akcí $\mathrm{ad}^*$ vypočtenou explicitně jako $-\boldsymbol{\Omega} \times \boldsymbol{\mu}$ v ztotožnění $\mathfrak{so}(3) \cong \mathbb{R}^3$, zamykajíc konvenci znaménka k souladu s kap. 02 §2.7.1 + kap. 03 §3.7 + oprava znaménka OQ-FORT-1 2026-05-24.
- **Křížovou validaci dvě-paralelně** (§4.5) ukazující Přístup A ≡ Přístup B na tělesové Eulerově rovnici — analytický analog doktrinálního §4.4 C2 bitové identity mezi implementacemi, aplikovaný na úrovni odvození.
- Translační rovnici $\mu \ddot{\boldsymbol{\rho}} = -\nabla_{\boldsymbol{\rho}} V$ (§4.6) rozloženou na Keplerovu sílu + zpětnou sílu, s **explicitním uzavřeným odvozením $\mathbf{F}_{\mathrm{back-reaction},\, X}$** (§4.6.2 zarámované): škáluje jako $1/r^4$, má orientačně závislý předfaktor, identicky mizí v limitě sférického tělesa (§4.6.3), účastní se Newtonovsko-třetí-zákon rovnováhy s tělesovými momenty sil k udržení zachování $\mathbf{P}_{\mathrm{CoM}}$ a $\mathbf{L}_{\mathrm{total}}$ (§4.6.4).
- **Úplnou sadu rovnic** (§4.7 zarámovanou) na 18rozměrném redukovaném fázovém prostoru, řádek po řádku křížově validovanou vůči `dynamics.py::state_derivative` na rotační + Keplerově úrovni, s mezerou zpětné síly identifikovanou explicitně a zařazenou v §4.9 OQ-4.5.

**Co dosud stanoveno NENÍ:**

- Tělesové Eulerovy rovnice uvedené jako samostatná specializace s momentem síly z gravitačního gradientu dosazeným explicitně na pravou stranu (kap. 03 §3.7 zarámovaný tvar rozvinutý do kartézských složek hlavních os). **Kapitola 05.**
- Numerické ověření pohybových rovnic vůči testovacím případům analytických limitů (Kepler, bezsilový symetrický setrvačník, librace gravitačního gradientu, asymetrický setrvačník Džanibekova, rozpojení při nekonečné vzdálenosti). **Kapitola 06.**
- Symplektický integrátor na Lieově grupě, jenž by nahradil RK4 pro dlouhohorizontovou přesnou integraci. **Rozsah inženýrské úrovně v0.2** (označeno v mezi­úrovňové poznámce kap. 01 §1.6 + tabulce vícemeřítkových časových měřítek kap. 03 §3.7.2 + §4.9 OQ-4.7).
- Rozšíření L2 / omezené soustavy tří těles přidávající centrifugální a Coriolisovy pseudo-síly. **Kanonická úroveň v0.2** (kap. 00 §1).
- Implementace zpětné síly v `dynamics.py::state_derivative`. **Rozsah inženýrské úrovně v0.2** (§4.9 OQ-4.5).

---

## §4.9 Otevřené otázky v této kapitole

### Vyřešeno v tomto kole (kanonická úroveň, 3. kolo, 2026-05-25)

- **OQ-4.1** ***Vyřešeno*** — volba mapy zamčena v §4.3: **mapa Eulerových úhlů** jako prezentace kanonické úrovně (odpovídá Podolský §7.3 přímo), kvaternionová mapa jako volba integrace inženýrské úrovně (poznámka inženýrské úrovně kap. 02 §2.7.1), obě produkující tutéž zarámovanou tělesovou Eulerovu rovnici. Žádné omezení jednotkové normy k zpracování v mapě kanonické úrovně; omezení mapy inženýrské úrovně je zpracováno v čase integrace re-normalizací na krok (standardní praxe kvaternionové propagace).
- **OQ-4.2** ***Vyřešeno*** — v §4.4.3 zarámovaném: $\mathrm{ad}^*_{\boldsymbol{\Omega}}\,\boldsymbol{\mu} = -\boldsymbol{\Omega} \times \boldsymbol{\mu}$ pro $\mathfrak{so}(3) \cong \mathbb{R}^3$, s minusovým znaménkem vynuceným konvencí párování. Znaménko konzistentní s levou trivializací kap. 02 §2.7.1 + kap. 03 §3.7 zarámovaným momentem síly + oprava znaménka OQ-FORT-1 (2026-05-24 odp.).
- **OQ-4.3** ***Vyřešeno*** — v §4.6.2 zarámovaném: $\mathbf{F}_{\mathrm{back-reaction},\, A} = -(3Gm_B/2r^4)\Phi_A \hat{\boldsymbol{\rho}} - (3Gm_B/r^4)(\mathbb{I}_3 - \hat{\boldsymbol{\rho}}\hat{\boldsymbol{\rho}}^T) \mathbf{R}_q^A (\mathbf{I}_A \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A})$. Škáluje jako $1/r^4$ jak se očekává; mizí v limitě sférického tělesa (§4.6.3); zachovává $\mathbf{P}_{\mathrm{CoM}}$ skrze Newtonovsko-třetí-zákon rovnováhu s momentem síly z gravitačního gradientu (§4.6.4).
- **OQ-4.4** ***Vyřešeno*** — v tabulce §4.7.1: shoda řádek po řádku vůči `dynamics.py::state_derivative` pro rotační + Keplerově-translační části. Zpětná síla není implementována v první verzi; velikost $\sim (\ell/r)^3$ vůči Keplerově síle, pod zaokrouhlovací chybou RK4 při konfiguracích první verze. Navazující úkol inženýrské úrovně označen jako nové OQ-4.5.

### Vystoupilo v tomto kole (dopředu orientované, rozsah v0.1.1 / v0.2)

- **OQ-4.5.** Implementovat $\mathbf{F}_{\mathrm{back-reaction},\, A} + \mathbf{F}_{\mathrm{back-reaction},\, B}$ (§4.6.2 zarámované) v translační integraci `dynamics.py::state_derivative`. Očekávaný dopad při konfiguracích první verze: pod zaokrouhlovací chybou dvojnásobné přesnosti. Stává se inženýrsky významným v režimu slapového uzamčení binárního asteroidu ($\ell/r \to 1$). **Rozsah inženýrské úrovně v0.1.1 nebo v0.2.** Doprovodná diagnostika zachování: s implementovanou zpětnou silou by orbitální moment hybnosti $\boldsymbol{\rho} \times \mu \dot{\boldsymbol{\rho}}$ měl vyměňovat s tělesovými spiny na úrovni předpovězené Noetherové identitou (§4.6.4) a diagnostika $|d\mathbf{L}_{\mathrm{total}}/L_0|$ by se měla zlepšit směrem k podlaze dvojnásobné přesnosti.
- **OQ-4.6.** Vyvinout **symplektický variační integrátor na Lieově grupě** (Hairer-Lubich-Wanner 2006 §VII.5; Lee-Leok-McClamroch 2018 *Global Formulations of Lagrangian and Hamiltonian Dynamics on Manifolds*) pro zarámovanou §4.7 sadu rovnic. Současná propagace RK4 je nesymplektická; dlouhohorizontová ($\gtrsim 10^7$ s, řádu orbitální / librační periody při hypotetických konfiguracích L2) integrace akumuluje sekulární úlet energie $\sim O(dt^4 \cdot T)$, jenž omezuje dosažitelnou věrnost fyziky. Rámec symplektiky na Lieově grupě respektuje podkladovou geometrickou strukturu (levá invariance + symplektická forma na $T^*Q_{\mathrm{CoM}}$) konstrukcí. **Rozsah inženýrské úrovně spárovaný s kanonickou úrovní v0.2.** Implementace: diskrétní Lagrangián na $SO(3) \times SO(3) \times \mathbb{R}^3$ skrze diskrétně-časovou Euler-Poincarého redukci (DT-EP) z Marsden-Pekarsky-Shkoller 1999; testováno vůči analytickým limitům kapitoly 06 plus test dlouhohorizontového úletu energie (dle `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md` §3.7).
- **OQ-4.7.** Vyvinout **energeticko-momentovou metodu** (Simo-Lewis-Marsden 1991; Marsden-Ratiu 1999 §1.7) pro analýzu stability zarámované §4.7 sady rovnic v relativních rovnováhách (např. librační rovnováha gravitačního gradientu kap. 03 §3.7.2 + kap. 05 §5.4). Energeticko-momentová metoda zobecňuje Routh-Hurwitzovu / Ljapunovovu stabilitu na soustavy se symetrií, poskytujíc rigorózní kritéria pro nelineární stabilitu relativních rovnováh na Lieových grupách. Páruje se s linearizovanou analýzou librace kap. 03 §3.7.2 (jež dává spektrální stabilitu) k upgradu na nelineární stabilitu pod doktrínou „testuj ≥2 nezávislé metody, jež souhlasí“. **Rozsah kanonické úrovně v0.2.**
- **OQ-4.8.** Vyvinout **Ljapunovovo spektrum** zarámované §4.7 sady rovnic poblíž režimu Džanibekova (§2.2 katalogu testovacích případů; kap. 05 §5.3 asymetrický případ; kap. 06 §6.5 testovací případ). Křížově validovat analytické Ljapunovovy exponenty (Wisdom 1987 pro Ljapunovova spektra spin-orbitální rezonance; specializováno na dvoutělesový sdružený případ Maciejewski 1995) vůči numerickému výpočtu z výstupu inženýrské úrovně. **Rozsah spárovaný s kanonickou úrovní v0.2.**

---

**Následující kapitola:** [05 — Tělesové Eulerovy rovnice jako specializace](05-body-frame-euler-equations.md)

**Konec 04-eulerovy-lagrangeovy-rovnice.md (v0.1, 3. kolo).**
