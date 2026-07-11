# Kapitola 6 — Analytické limity a validační brány

> **Předpokládané čtení.** Kapitoly 01–05. Tato kapitola pracuje se zarámovanou celkovou sadou rovnic kapitoly 04 §4.7 + specializací na hlavní osy z kapitoly 05, validujíc každou část vůči konfiguracím analytických limitů zkatalogizovaným v `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md`.
>
> **Stav.** v0.1 (kanonická úroveň, 3. kolo, 2026-05-25). Lešení + obsah pro testovací případy §6.1–§6.7; §6.8 stanovení + §6.9 seznam OQ pro práci 4. kola+.

> **Poznámka ke způsobu a předpokladům odvození.** Tato CS kanonická kapitola (v0.1) je odvozena z anglické kanonické kapitoly `docs/canonical/en/06-analytical-limits-and-validation-gates.md` (v0.1, 3. kolo). EN verze zůstává autoritativním pramenem při budoucích aktualizacích.

---

## §6.1 Co tato kapitola dělá

Předchozích pět kapitol vybudovalo pohybové rovnice kanonické úrovně z prvních principů: konfigurační varieta + variační princip (kap. 01) → kinetická energie jako levo-invariantní Riemannova metrika (kap. 02) → vzájemný gravitační potenciál skrze multipolový rozvoj (kap. 03) → sestavení Lagrangiánu + Euler-Poincarého redukce (kap. 04) → Eulerovy rovnice hlavních os (kap. 05). Výsledkem je zarámovaná sada rovnic kap. 04 §4.7 na 18rozměrném redukovaném fázovém prostoru $T^* Q_{\mathrm{CoM}}$, s rotační dynamikou specializovanou na tvar hlavních os v kap. 05 §5.4 zarámovaném.

Odvození kanonické úrovně, jež je vnitřně konzistentní, avšak nefalzifikovatelné, by si nezasloužilo své místo pod lege-artis disciplínou. Kapitola 06 přivádí odvozené rovnice do kontaktu s **falzifikovatelnými testovacími případy analytických limitů**: konfiguracemi soustavy, pro něž pohybové rovnice připouštějí uzavřená analytická řešení, vůči nimž lze jakoukoli numerickou implementaci kvantitativně srovnat. Disciplínou je kritérium doktríny §4.4 C1 (učebnicově-doslovné pokrytí testy ≥3 případy na publikovaný anchor), přizpůsobené doméně analytické mechaniky.

Podle konvence zamčené v `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md` §1 se testovací případy dělí do dvou tříd:

- **Třída A — Přesná uzavřená analytická řešení.** Konfigurace, kde se pohybové rovnice redukují na soustavu, jejíž řešení je známé v uzavřené formě. Analytické řešení slouží jako **orákulum**: v libovolném čase $t$ je analytická hodnota pravdou, vzdálenost numerické hodnoty od ní je měřitelná chyba. Třída A odpovídá na otázku *„řeší integrátor rovnice správně?“*
- **Třída B — Konfigurace okrajových případů.** Konfigurace, kde se dynamika buď redukuje na jednodušší problém v jistém asymptotickém režimu, nebo se blíží singulární hranici oboru platnosti modelu. Testovací případy třídy B ověřují, že integrátor zvládá přechody elegantně a že jsou respektovány uvedené předpoklady modelu. Třída B odpovídá na otázku *„řeším správné rovnice?“*

Obě třídy jsou nutné, neboť odpovídají na různé otázky. Bez třídy B může být model numericky přesný, avšak fyzikálně chybný. Bez třídy A platí opak.

Každý testovací případ níže je dodán s: matematickým nastavením, odvozením redukce analytického limitu nebo očekávaného chování, uzavřeným řešením (třída A) nebo očekávaným chováním (třída B), CS-jazykovým ukotvením Podolský 2018 tam, kde je dostupné, očekávanou numerickou tolerancí a ukazatelem na umístění testovacího souboru inženýrské úrovně / Fortran fáze 1. Katalog v0.1 je strukturován kolem **šesti základních testovacích případů** §6.2–§6.6 plus procedurální §6.7 brány mezi implementacemi. v0.2 rozšíří o rodinu omezené soustavy tří těles u L2.

## §6.2 Testovací případ Keplerovy limity (třída A)

### §6.2.1 Matematické nastavení

Nastav oba tenzory setrvačnosti na nulu: $\mathbf{I}_A \to \mathbf{0}, \mathbf{I}_B \to \mathbf{0}$. Podle kap. 03 §3.4.3 zarámovaného obě slapové části $V_{\mathrm{tidal},A}, V_{\mathrm{tidal},B}$ identicky mizí (každá je úměrná stopě setrvačnosti + bilineární projekci setrvačnosti). Plný vzájemný potenciál se zhroutí na monopólově-monopólový Keplerův člen:

$$
V \;\xrightarrow{\mathbf{I}_{A,B} \to \mathbf{0}}\; -\frac{G m_A m_B}{|\boldsymbol{\rho}|}.
$$

Podle numerické kontroly kap. 04 §4.7.2 se rotační rovnice rozpojí (oba $\boldsymbol{\tau}_X^{\mathrm{body}} = \mathbf{0}$ a zpětná síla $\mathbf{F}_{\mathrm{back-reaction},\, X} = \mathbf{0}$ mizí, neboť obě škálují s $\mathbf{I}_X$). Zarámovaná §4.7 sada rovnic se redukuje na:

$$
\mu \ddot{\boldsymbol{\rho}} \;=\; -\frac{Gm_A m_B}{|\boldsymbol{\rho}|^3}\, \boldsymbol{\rho}, \qquad \dot{\boldsymbol{\Omega}}_X = \mathbf{0}.
$$

Translační rovnice je standardním Keplerovým problémem s redukovanou hmotností; rotační rovnice se rozpojí na $\boldsymbol{\Omega}_X = \mathrm{const}$.

### §6.2.2 Analytické řešení

V těžišťové soustavě relativní souřadnice $\boldsymbol{\rho}$ vykresluje kuželosečku (elipsu, parabolu, hyperbolu) určenou energií $E$ a momentem hybnosti $|\mathbf{L}|$ relativního pohybu. Pro vázanou (eliptickou) orbitu s hlavní poloosou $a$ je orbitální perioda:

$$
T_{\mathrm{orbit}} \;=\; 2\pi \sqrt{\frac{a^3}{G(m_A + m_B)}},
$$

(třetí Keplerův zákon). Orbitální energie je $E = -G m_A m_B / (2a)$, excentricita je $e = \sqrt{1 + 2 E |\mathbf{L}|^2 / (\mu (Gm_Am_B)^2)}$ a Laplaceův-Runge-Lenzův vektor $\mathbf{A} = \mathbf{p}_{\mathrm{rel}} \times \mathbf{L} - \mu Gm_Am_B \hat{\boldsymbol{\rho}}$ je přesně zachován (jeho zachování odráží skrytou $SO(4)$ symetrii Keplerova problému; viz Landau-Lifšic, sv. I §15).

### §6.2.3 Podolského ukotvení (CS-jazykové)

Keplerův problém je pojednán v klasických česko-jazykových učebnicích analytické mechaniky, avšak není specificky ukotven v Podolský 2018 Kapitola 6–8 (jež se zaměřují na kinematiku + dynamiku + aplikace tuhého tělesa). Pro CS čtenáře je standardním precedentem Landau-Lifšic, sv. I §13–§15 (CS překlad *Mechaniky* je široce dostupný v českých univerzitních knihovnách).

### §6.2.4 Numerická tolerance + testovací soubor

| Veličina | Tolerance | Testovací soubor |
|---|---|---|
| Orbitální perioda vs $2\pi\sqrt{a^3/(GM)}$ | $< 10^{-6}$ relativně | `tests/test_dynamics.py::TestKeplerLimit` (plánováno) |
| Stabilita excentricity přes $T_{\mathrm{orbit}}$ | $< 10^{-6}$ | (totéž) |
| Velikost LRL vektoru přes $T_{\mathrm{orbit}}$ | $< 10^{-8}$ | (totéž) |
| Zachování $E$ a $|\mathbf{L}|$ | $< 10^{-10}$ | (totéž) |

Status implementace: plánováno (dle `_specs/...` §2.3). Cesta inženýrské úrovně: nastav `geometries.make_zero_inertia_body()` (syntetický konstruktor, zařazeno) k vyprodukování zástupců s nulovou setrvačností; vyvolej `dynamics.py::integrate` s jinak-standardní konfigurací; srovnej vůči referenci Keplerovy rovnice `scipy.integrate.solve_ivp` při stejném $dt$.

## §6.3 Testovací případ volného symetrického setrvačníku (třída A)

### §6.3.1 Matematické nastavení

Jediné tuhé těleso s axisymetrickým tenzorem setrvačnosti $\mathbf{I} = \mathrm{diag}(I_\perp, I_\perp, I_\parallel)$; žádný vnější moment síly ($\boldsymbol{\tau}^{\mathrm{body}} = \mathbf{0}$). Zarámované Eulerovy rovnice kap. 05 §5.2.1 se specializují (dle kap. 05 §5.3.2) na:

$$
I_\perp \dot\Omega_x = (I_\perp - I_\parallel)\,\Omega_y \Omega_z, \quad
I_\perp \dot\Omega_y = -(I_\perp - I_\parallel)\,\Omega_x \Omega_z, \quad
I_\parallel \dot\Omega_z = 0.
$$

Bezsilový případ je v našem dvoutělesovém nastavení dosažitelný (a) nastavením hmotnosti druhého tělesa na nulu, takže není žádný Keplerův partner a žádné pole gravitačního gradientu, nebo (b) umístěním druhého tělesa do nekonečné vzdálenosti, takže slapová vazba mizí (kap. 03 §3.8.4 / §6.6 níže). Implementačně nejčistší cestou inženýrské úrovně je (a) s jednotělesovým integračním vstupem, což je to, co testbed fáze 1 D8 `backends/fortran/tests/test_symmetric_top.f90` implementuje.

### §6.3.2 Analytické řešení

$\Omega_z(t) = \Omega_z(0) = \mathrm{const}$ (zachováno). Příčné složky splňují vázanou harmonickou oscilaci:

$$
\Omega_x(t) = A \cos(\lambda t + \phi_0), \qquad
\Omega_y(t) = A \sin(\lambda t + \phi_0), \qquad
\Omega_z(t) = \Omega_z^0,
$$

s **Eulerovou precesní úhlovou frekvencí** (často zvanou „rychlost precese polovičního úhlu tělesového kuželu“):

$$
\lambda \;=\; \frac{I_\parallel - I_\perp}{I_\perp}\, \Omega_z^0.
$$

Vektor tělesové úhlové rychlosti $\boldsymbol{\Omega}(t)$ vykresluje **tělesový kužel** polovičního úhlu $\arctan(A/\Omega_z^0)$ kolem osy symetrie $\hat{\mathbf{e}}_3$ rychlostí $\lambda$. Zachovávající se veličiny: kinetická energie $T = \tfrac{1}{2}(I_\perp A^2 + I_\parallel \Omega_z^{0\,2})$ a velikost momentu hybnosti $|\mathbf{L}|^2 = I_\perp^2 A^2 + I_\parallel^2 \Omega_z^{0\,2}$ obě přesně zachovány.

### §6.3.3 Podolského ukotvení (CS-jazykové)

**Přímé ukotvení.** Podolský 2018 §8.1 *Volný symetrický setrvačník* odvozuje toto přesné řešení z jeho zarámované Rovn. (7.15) bezsilové + axisymetrické setrvačnosti. Jeho Rovn. (8.1) je analogem našeho řešení $(\Omega_x, \Omega_y, \Omega_z)$ výše a jeho Rovn. (8.2) dává geometrickou interpretaci tělesového kuželu. Fortranový testbed fáze 1 D8 `backends/fortran/tests/test_symmetric_top.f90` užívá Podolský Rovn. (8.1)–(8.2) jako CS-jazykové analytické orákulum ve srovnání.

### §6.3.4 Numerická tolerance + testovací soubor

| Veličina | Tolerance | Testovací soubor |
|---|---|---|
| Úlet $\Omega_z$ přes $T = 600$ s | $< 10^{-12}$ relativně | `backends/fortran/tests/test_symmetric_top.f90` (fáze 1 D8) + `tests/test_dynamics.py::TestRK4OnTorqueFreeSymmetricTop` (inženýrská úroveň, implementováno v první verzi) |
| Stabilita obálky $A = |\Omega_x(t) + i\Omega_y(t)|$ | $< 10^{-12}$ relativně | (totéž) |
| Chyba obnovení precesní frekvence $\lambda$ | $< 0{,}2$ FFT binů (omezeno vzorkováním, nikoli integrací) | (totéž) |
| Zachování energie / momentu hybnosti | $|dE/E_0| < 10^{-11}$, $|d\mathbf{L}/L_0| < 10^{-9}$ | (totéž) |

**Status: implementováno v první verzi** na inženýrské úrovni (brána zachování `TestRK4OnTorqueFreeSymmetricTop` nyní prochází při $6{,}2 \times 10^{-12}$ pro energii a $2{,}1 \times 10^{-10}$ pro moment hybnosti; kap. 00 §3 cituje tyto jako hlavní diagnostiky první verze). Fortranový křížový port fáze 1 D8 je lege-artis křížovou kontrolou bitové identity.

## §6.4 Testovací případ librace gravitačního gradientu (třída A linearizovaná + třída A nelineární přes redukci Lagrangeova setrvačníku)

### §6.4.1 Matematické nastavení

Axisymetrické těleso $A$ ($I_{xx,A} = I_{yy,A} \equiv I_\perp$, $I_{zz,A} \equiv I_\parallel$) na kruhové orbitě kolem bodové hmoty tělesa $B$ ($\mathbf{I}_B = \mathbf{0}$), s osou symetrie tělesa $\hat{\mathbf{e}}_3^A$ nominálně srovnanou s $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$ (radiálně směřující rovnováha). Zarámované dosazené Eulerovy rovnice kap. 05 §5.4 se specializují s $\boldsymbol{\tau}_A^{\mathrm{body}}$ z gravitačního gradientu + axisymetrickým $\mathbf{I}_A$. Orbitální rychlost je $n = \sqrt{G(m_A + m_B)/|\boldsymbol{\rho}|^3}$.

Perturbuj orientaci o malé úhly $(\theta, \phi)$ kolem rovnováhy (perturbace klopení–klonění, s $\theta$ rotací kolem tělesové $\hat{\mathbf{e}}_1$ a $\phi$ kolem tělesové $\hat{\mathbf{e}}_2$). Do lineárního řádu v $(\theta, \phi)$ se tělesový směr vzdálenosti stává $\hat{\boldsymbol{\rho}}_{\mathrm{body},A} \approx (-\phi, -\theta, 1)$ a dosazené rovnice kap. 05 §5.4 dávají linearizované librační rovnice (odvození kap. 05 §5.4.1 rozšířené na oba příčné módy; zarámovaná orientační Jacobiho matice kap. 03 §3.7.2 poskytuje ekvivalentní maticovou formu).

### §6.4.2 Uzavřené řešení — linearizovaný režim

Linearizovaná librace je dvourozměrná, s frekvencemi vlastních módů určenými poměrem setrvačnosti $(I_\parallel - I_\perp)/I_\perp$ a orbitální rychlostí $n$. Pro rovinný klopivý mód kolem radiálně směřující rovnováhy (odvození kap. 05 §5.4.1; křížově kontrolováno vůči orientační Jacobiho matici kap. 03 §3.7.2):

$$
\omega_{\mathrm{lib}}^2 \;=\; \frac{3 G m_B (I_\perp - I_\parallel)}{I_\perp\, |\boldsymbol{\rho}|^3}
\;=\; 3 n^2\, \frac{I_\perp - I_\parallel}{I_\perp},
$$

s užitím $n^2 = G m_B / |\boldsymbol{\rho}|^3$ pro orbitální rychlost v limitě bodové hmoty primáru (kap. 06 §6.2 Keplerova limita). Pro protáhlé ($I_\parallel < I_\perp$): $\omega_{\mathrm{lib}}^2 > 0$ → stabilní librace. Pro zploštělé ($I_\parallel > I_\perp$): $\omega_{\mathrm{lib}}^2 < 0$ → exponenciální nestabilita; postoj s osou symetrie podél radiály je nestabilní rovnováhou a těleso relaxuje k postoji s osou symetrie kolmou k radiále. Plné pojednání trojmódové vazby klopení-klonění-zatáčení je v Hughes 2004 §8.4 Rovn. (8.40)–(8.46).

### §6.4.3 Uzavřené řešení — plně nelineární přes efektivní potenciál Lagrangeova setrvačníku (řeší OQ-5.6 + OQ-6.4)

Librace o konečné amplitudě kolem radiálně směřující rovnováhy je analogem klasického problému **Lagrangeova setrvačníku** (těžký symetrický setrvačník s pevným bodem v homogenní gravitaci) v dynamice postoje kosmické lodě. Spojení je strukturální: axisymetrické těleso v polohově závislém gravitačním poli s osou symetrie nakloněnou o úhel $\vartheta$ od směru pole splňuje Lagrangián téhož algebraického tvaru bez ohledu na to, zda je „pole“ homogenním gravitačním zrychlením (Podolský §8.3) nebo polohově závislým momentem síly z gravitačního gradientu (naše nastavení §6.4). Odvození Podolský 2018 §8.3 se přenáší doslova se substitucí $g \to 3 G m_B / |\boldsymbol{\rho}|^3$ jako parametrem efektivního zrychlení a **metoda efektivního potenciálu** Lagrangeova setrvačníku poskytuje uzavřené řešení o konečné amplitudě.

**Nastavení.** Přijmi intrinzickou-ZXZ mapu Eulerových úhlů kap. 04 §4.3.1 / Podolský §6.4 Rovn. (6.14)–(6.16) s osou symetrie tělesa-$A$ podél tělesové $\hat{\mathbf{e}}_3$. V radiálně směřující rovnováze jsou Eulerovy úhly $(\varphi, \vartheta, \psi) = (\varphi_0(t), 0, \psi_0(t))$, s $\varphi$ kódujícím orbitální pohyb a $\psi$ spin kolem osy symetrie. Perturbace mimo rovnováhu dává $\vartheta \neq 0$ — **nutační úhel** mezi tělesovou $\hat{\mathbf{e}}_3$ a $\hat{\boldsymbol{\rho}}$.

Lagrangián pro axisymetrické těleso v poli gravitačního gradientu (kap. 04 §4.2 zarámovaný, specializovaný na axisymetrickou setrvačnost a mapu Eulerových úhlů kap. 04 §4.3) zní:

$$
L \;=\; \tfrac{1}{2} I_\perp (\dot\vartheta^2 + \dot\varphi^2 \sin^2\vartheta)
\;+\; \tfrac{1}{2} I_\parallel (\dot\psi + \dot\varphi \cos\vartheta)^2
\;-\; V_{\mathrm{tidal},\, A}(\vartheta).
$$

V režimu, kde je orbitální poloměr velký vůči rozměru tělesa a orbitální perioda mnohem delší než perioda spinu — obojí dobře splněno v hypotetických konfiguracích podobných L2 (tabulka vícemeřítkových časových měřítek kap. 03 §3.7.2) — se slapový potenciál stává efektivně orientačně závislým jen skrze $\vartheta$:

$$
V_{\mathrm{tidal},\, A}(\vartheta)
\;=\;
-\frac{G m_B}{2 |\boldsymbol{\rho}|^3}\,
\big[\, \mathrm{tr}(\mathbf{I}_A) - 3 (I_\perp \sin^2\vartheta + I_\parallel \cos^2\vartheta)\,\big]
\;=\;
-\frac{G m_B}{2 |\boldsymbol{\rho}|^3}\,
\big[\, (I_\perp + I_\parallel) - 3 I_\perp + 3 (I_\perp - I_\parallel) \cos^2\vartheta \,\big],
$$

až na konstantu nezávislou na $\vartheta$ absorbovanou do energetické reference. Orientačně závislá část je **$\propto \cos^2\vartheta$**, což je přesně algebraický tvar Podolského gravitačního potenciálu v nastavení Lagrangeova setrvačníku (jeho Rovn. 8.7–8.9 se substitucí poznamenanou výše). Lagrangián má tedy **cyklické souřadnice** $\varphi$ a $\psi$ — ani jedna se v $L$ neobjevuje explicitně — dávaje dva první integrály skrze Noetherové větu (kap. 02 §2.6.1 + rámec zachování kap. 04 §4.6.4):

$$
L_\varphi \;\equiv\; \frac{\partial L}{\partial \dot\varphi}
\;=\; I_\perp \dot\varphi \sin^2\vartheta + I_\parallel (\dot\psi + \dot\varphi \cos\vartheta)\cos\vartheta
\;=\; \mathrm{const},
$$

$$
L_\psi \;\equiv\; \frac{\partial L}{\partial \dot\psi}
\;=\; I_\parallel (\dot\psi + \dot\varphi \cos\vartheta)
\;=\; \mathrm{const}.
$$

Toto jsou **azimutální složka momentu hybnosti** (kolem $\hat{\boldsymbol{\rho}}$) a **tělesový spinový moment hybnosti** (kolem tělesové $\hat{\mathbf{e}}_3$), po řadě — přesně Podolského Rovn. (8.10)–(8.11). Se zachovanými $L_\psi$ + $L_\varphi$ se pohybová rovnice pro zbývající stupeň volnosti $\vartheta$ redukuje kvadraturou: dosaď $\dot\psi$ a $\dot\varphi$ z prvních integrálů zpět do celkové mechanické energie $E = T + V$ (třetí zachovávající se veličiny) a vyřeš $\dot\vartheta^2$ jako funkci samotného $\vartheta$.

**Efektivní potenciál.** Proveď substituci. Z $L_\psi = I_\parallel(\dot\psi + \dot\varphi\cos\vartheta)$ máme $\dot\psi + \dot\varphi\cos\vartheta = L_\psi/I_\parallel$, a tak spinová kinetická energie $\tfrac{1}{2}I_\parallel(\dot\psi + \dot\varphi\cos\vartheta)^2 = L_\psi^2/(2I_\parallel) = \mathrm{const}$ — absorbována do energetické reference. Z $L_\varphi = I_\perp\dot\varphi\sin^2\vartheta + L_\psi\cos\vartheta$ řešíme $\dot\varphi = (L_\varphi - L_\psi\cos\vartheta)/(I_\perp\sin^2\vartheta)$. Dosazením do rotační kinetické energie:

$$
T_{\mathrm{rot,\, effective}}
\;=\;
\tfrac{1}{2} I_\perp \dot\vartheta^2
\;+\;
\frac{(L_\varphi - L_\psi \cos\vartheta)^2}{2 I_\perp \sin^2\vartheta}
\;+\; \mathrm{const}.
$$

První člen je kinetickou energií stupně volnosti $\vartheta$; druhý je **centrifugálně-bariéře podobným** příspěvkem z eliminovaných cyklických souřadnic. Definujíc **efektivní potenciál** $V_{\mathrm{ef}}(\vartheta)$, jenž kombinuje tento centrifugální člen s orientačně závislou částí gravitačního gradientu:

$$
\boxed{
V_{\mathrm{ef}}(\vartheta)
\;=\;
\frac{(L_\varphi - L_\psi \cos\vartheta)^2}{2 I_\perp \sin^2\vartheta}
\;-\;
\frac{3 G m_B (I_\perp - I_\parallel)}{2 |\boldsymbol{\rho}|^3}\, \cos^2\vartheta.
}
$$

Plná jednorozměrná pohybová rovnice v $\vartheta$ je pak:

$$
\tfrac{1}{2} I_\perp \dot\vartheta^2 \;+\; V_{\mathrm{ef}}(\vartheta) \;=\; E_{\mathrm{ef}} \;=\; \mathrm{const}.
$$

Toto je **redukce Lagrangeova setrvačníku**: rotační problém o dvou stupních volnosti (po vyintegrování spinu) redukovaný na jediný jednorozměrný pohyb $\vartheta$ v efektivním potenciálu. Podolského Rovn. (8.12)–(8.15) dávají přesně tuto strukturu s naší substitucí.

**Uzavřené řešení kvadraturou.** Trajektorie $\vartheta(t)$ se získá separací proměnných:

$$
\boxed{
t \;-\; t_0 \;=\;
\int_{\vartheta_0}^{\vartheta(t)}\, \frac{d\vartheta'}{\sqrt{2(E_{\mathrm{ef}} - V_{\mathrm{ef}}(\vartheta'))/I_\perp}}.
}
$$

Integrál je obecně **eliptický** (kubická algebraická rovnice $E_{\mathrm{ef}} = V_{\mathrm{ef}}(\vartheta)$ pro body obratu má tři kořeny $\vartheta_1, \vartheta_2, \vartheta_3$; integrál se redukuje na eliptické integrály prvního a třetího druhu, vyjádřené v pojmech Jacobiho eliptických funkcí nebo ekvivalentně Weierstrassovy funkce $\wp$). Při malé amplitudě — $|\vartheta| \ll 1$ — se eliptický integrál redukuje na linearizovanou libraci z §6.4.2: $\vartheta(t) = \vartheta_{\mathrm{amp}} \cos(\omega_{\mathrm{lib}} t + \phi_0)$ s linearizovaným vzorcem frekvence obnoveným přesně. Při konečné amplitudě se perioda **prodlužuje** s amplitudou (škálování perioda-amplituda: $T_{\mathrm{lib}}(\vartheta_{\mathrm{amp}})/T_{\mathrm{lib}}(0) = (2/\pi) K(k)$, kde $K$ je úplný eliptický integrál prvního druhu a $k$ je eliptický modul určený $\vartheta_{\mathrm{amp}}$, spinem $L_\psi$ a poměrem setrvačnosti).

**Tři precesní režimy (Podolský §8.3 (a)/(b)/(c)).** V závislosti na poměru $L_\psi / L_\varphi$ a energii $E_{\mathrm{ef}}$ jáma efektivního potenciálu připouští tři kvalitativně odlišné trajektorie na jednotkové kouli (hrot tělesové $\hat{\mathbf{e}}_3$):

- **Režim (a) — čistá precese:** $\vartheta = \mathrm{const}$ (tělesová $\hat{\mathbf{e}}_3$ vykresluje kružnici na kouli). Podmínkou je $dV_{\mathrm{ef}}/d\vartheta = 0$ ve zvolené hodnotě $\vartheta$; fyzikálně spinem indukovaná centrifugální bariéra přesně vyvažuje moment síly z gravitačního gradientu.
- **Režim (b) — nutace mezi dvěma zeměpisnými šířkami $\vartheta_1 < \vartheta_2$:** standardní librační režim; $\vartheta$ osciluje mezi body obratu $\vartheta_1$ a $\vartheta_2$, kde $V_{\mathrm{ef}}(\vartheta_i) = E_{\mathrm{ef}}$. Tělesová $\hat{\mathbf{e}}_3$ vykresluje vlnitou dráhu mezi dvěma rovnoběžkami na jednotkové kouli; hroty se objevují v bodech obratu, když je tam $\dot\varphi = 0$.
- **Režim (c) — smyčkující nutace:** je-li $L_\psi < L_\varphi/(1+|\cos\vartheta_{\max}|)$ nebo podobná prahová podmínka, azimutální úhlová rychlost $\dot\varphi$ mění znaménko v rámci nutačního cyklu; trajektorie ukazuje smyčky na kouli.

Numerické škálování amplitudy librační periody $T_{\mathrm{lib}}(\vartheta_{\mathrm{amp}})$ — **korekce o konečné amplitudě** k linearizovanému vzorci §6.4.2 — je klíčovou experimentálně falzifikovatelnou predikcí uzavřeného řešení, s explicitní závislostí na eliptickém modulu.

### §6.4.4 Podolského ukotvení (CS-jazykové) — varianta Lagrangeova setrvačníku

**Přímé ukotvení.** Podolský 2018 §8.3 *Těžký symetrický setrvačník s pevným bodem* (Lagrangeův setrvačník) poskytuje CS-jazykové analytické orákulum. Konkrétně:

- Podolský Rovn. (8.10)–(8.11) jsou naše první integrály cyklických souřadnic $L_\psi, L_\varphi$ — identická algebraická struktura se substitucí $g \to 3Gm_B/|\boldsymbol{\rho}|^3$.
- Podolský zarámovaná Rovn. (8.15) je náš §6.4.3 zarámovaný efektivní potenciál $V_{\mathrm{ef}}(\vartheta)$ — člen centrifugální bariéry identický; gravitační člen substituuje homogenní-gravitaci $V = mgl\cos\vartheta$ naším $-3Gm_B(I_\perp-I_\parallel)\cos^2\vartheta/(2|\boldsymbol{\rho}|^3)$. Funkční tvar se liší (jeho $\cos\vartheta$ vs naše $\cos^2\vartheta$), neboť Podolského nastavení má *těžiště* tělesa posunuté od pevného bodu (dávaje kosinovou závislost na úhlu náklonu od gravitace), zatímco naše má těleso *nakloněné* od radiálního směru v polohově závislém gravitačním poli (dávaje kosinově-kvadratickou závislost). Kvalitativní fenomenologie — omezená nutace, tři režimy (a)/(b)/(c) — je identická.
- Podolský Rovn. (8.13) s případy (a)/(b)/(c) mapuje přímo na naše tři režimy z §6.4.3.

Testbed fáze 2 `backends/fortran/tests/test_libration.f90` (zařazeno pro fázi 2 fortranového lege-artis portu) bude křížově validovat vůči zarámovanému §6.4.3 vzorci $V_{\mathrm{ef}}$ numericky: integruj zarámovanou sadu rovnic kap. 04 §4.7 s axisymetrickým tělesem v poli gravitačního gradientu, extrahuj librační periodu při více amplitudách a srovnej vůči analytické periodě $T_{\mathrm{lib}}(\vartheta_{\mathrm{amp}}) = (2/\pi) K(k(\vartheta_{\mathrm{amp}})) \cdot T_{\mathrm{lib}}^{\mathrm{linear}}$. **(Řeší OQ-5.6 + OQ-6.4 společně.)**

### §6.4.5 Numerická tolerance + testovací soubor

| Veličina | Tolerance | Testovací soubor |
|---|---|---|
| Linearizovaná librační perioda vs $2\pi/\omega_{\mathrm{lib}}$ | $< 10^{-4}$ relativně přes 10 orbitálních period | `backends/fortran/tests/test_libration.f90` (fáze 2; zařazeno) + `tests/test_integration.py::TestGravityGradientLibration` (plánovaná inženýrská úroveň) |
| Librační perioda o konečné amplitudě vs $\boxed{T_{\mathrm{lib}}(\vartheta_{\mathrm{amp}}) = (2/\pi) K(k) \cdot T_{\mathrm{lib}}^{\mathrm{linear}}}$ (§6.4.3) | $< 10^{-3}$ relativně při amplitudách do $\vartheta_{\mathrm{amp}} \lesssim \pi/4$; systematicky se zhoršuje při větších amplitudách, jak se rozvoj eliptického modulu láme | (totéž; sken konečné amplitudy fáze 2) |
| Stabilita amplitudy (žádné nefyzikální tlumení nebo růst) | $< 10^{-6}$ relativně přes 10 orbitálních period | (totéž) |
| Klasifikace stability (protáhlé-stabilní vs zploštělé-nestabilní) | Znaménkově-správná z vlastní hodnoty orientační Jacobiho matice D5 (regresní pojistka) | `backends/fortran/tests/test_orientation_jacobian.f90` (fáze 1 D5) |
| Klasifikace tří režimů (a)/(b)/(c) | Trajektorie na jednotkové kouli hrotu tělesové $\hat{\mathbf{e}}_3$ odpovídá Podolský §8.3 fenomenologii při odpovídajících poměrech $L_\psi / L_\varphi$ | Doplňkové fáze 2 |

S uzavřeným řešením z §6.4.3 je celý testbed třídy A (linearizovaný režim i režim o konečné amplitudě nyní mají analytická orákula). Předchozí označení „třída B pouze-numerická“ pro režim o konečné amplitudě je **retirováno** — redukce Lagrangeova setrvačníku poskytuje chybějící analytický obsah.

## §6.5 Testovací případ asymetrického setrvačníku Džanibekova (třída A linearizovaná + třída A nelineární přes Jacobiho eliptické funkce)

### §6.5.1 Matematické nastavení

Jediné tuhé těleso se třemi odlišnými hlavními momenty $I_{xx} < I_{yy} < I_{zz}$; žádný vnější moment síly. Nastav počáteční úhlovou rychlost nominálně kolem prostřední osy $\hat{\mathbf{e}}_2$: $\boldsymbol{\Omega}(0) = (0, \Omega_2^0, 0)$ + malá perturbace $(\delta\Omega_x, \delta\Omega_y, \delta\Omega_z)$ s $|\delta\boldsymbol{\Omega}| \ll |\Omega_2^0|$. Zarámované Eulerovy rovnice kap. 05 §5.2.1 linearizované kolem tohoto počátečního stavu dávají (odvození kap. 05 §5.3.3):

$$
I_{xx} \delta\dot\Omega_x = (I_{yy} - I_{zz})\, \delta\Omega_z\, \Omega_2^0, \qquad
I_{zz} \delta\dot\Omega_z = (I_{xx} - I_{yy})\, \delta\Omega_x\, \Omega_2^0.
$$

(Rovnice $\delta\Omega_y$ je v lineárním řádu rozpojena; perturbace zpočátku zachovává $\delta\Omega_y \approx 0$ až do nelineárního režimu.) Dosazením jedné rovnice do druhé:

$$
\delta\ddot\Omega_x \;=\; -\frac{(I_{yy} - I_{xx})(I_{zz} - I_{yy})}{I_{xx} I_{zz}}\, (\Omega_2^0)^2\, \delta\Omega_x.
$$

Koeficient vpravo je **záporný**, neboť $(I_{yy} - I_{xx}) > 0$ a $(I_{zz} - I_{yy}) > 0$ pro předpokládané uspořádání, takže závorkovaný čitatel je kladný a vedoucí znaménko činí rovnici $\delta\ddot\Omega_x = +\sigma^2 \delta\Omega_x$ s $\sigma^2 > 0$. Řešení je exponenciálně narůstající: $\delta\Omega_x \sim e^{\sigma t}$.

### §6.5.2 Uzavřená linearizovaná rychlost růstu

Linearizovaná rychlost růstu je:

$$
\sigma \;=\; \Omega_2^0\, \sqrt{\frac{(I_{yy} - I_{xx})(I_{zz} - I_{yy})}{I_{xx} I_{zz}}}.
$$

Toto charakterizuje raný exponenciální růst $|\delta\Omega_x(t)| \sim e^{\sigma t}$ předtím, než nelinearity nasytí perturbaci. Křížová validace s učebnicovou klasifikací teorému tenisové rakety (kap. 05 §5.2.2) je přesná.

### §6.5.3 Uzavřené řešení o konečné amplitudě přes Jacobiho eliptické funkce (řeší OQ-5.5)

Pro dynamiku o **konečné amplitudě** — plný vzor přeskoků Džanibekova s neinfinitezimální amplitudou perturbace — je bezsilový problém asymetrického setrvačníku **integrabilní kvadraturou**, neboť zarámované Eulerovy rovnice kap. 05 §5.2.1 připouštějí *dvě* zachovávající se veličiny (kap. 05 §5.5.1): tělesovou kinetickou energii $T$ a čtverec momentu hybnosti $|\mathbf{L}|^2$. Spolu s tělesovými Eulerovými rovnicemi tato dvě omezení redukují trojrozměrnou soustavu na jeden efektivní stupeň volnosti, jehož vývoj je dán v uzavřené formě **Jacobiho eliptickými funkcemi**.

**Redukce zákony zachování.** Označ $L_i \equiv I_i \Omega_i$ (tělesové složky momentu hybnosti; nezaměňovat s inerciálním celkovým momentem hybnosti). Dvě zachovávající se veličiny jsou:

$$
2T \;=\; L_x^2/I_1 + L_y^2/I_2 + L_z^2/I_3 \;=\; \mathrm{const}, \qquad
|\mathbf{L}|^2 \;=\; L_x^2 + L_y^2 + L_z^2 \;=\; \mathrm{const}.
$$

V tělesovém $\mathbf{L}$-prostoru leží trajektorie na průniku:
- **Koule momentu hybnosti** $L_x^2 + L_y^2 + L_z^2 = |\mathbf{L}|^2$ (koule o poloměru $|\mathbf{L}|$).
- **Elipsoidu kinetické energie** $L_x^2/I_1 + L_y^2/I_2 + L_z^2/I_3 = 2T$ (elipsoidu s poloosami $\sqrt{2 T I_i}$ podél hlavních os).

Genericky (pro $|\mathbf{L}|^2 \neq 2 T I_i$ pro žádné $i$) je průnikem hladká uzavřená křivka. Trajektorie prochází touto křivkou periodicky — **polhodie** v tělesovém $\mathbf{L}$-prostoru (názvosloví Goldstein 2002 §5.6).

**Parametrizace eliptickými funkcemi.** Předpokládej $I_1 < I_2 < I_3$ a $|\mathbf{L}|^2 > 2 T I_2$ (větev, jejíž polhodie obkružuje hlavní osu **největšího** momentu $\hat{\mathbf{e}}_3$ — tuto větev realizují testovací konfigurace §6.5.5, např. $\mathbf{I} = \mathrm{diag}(100,200,400)$, $\boldsymbol{\Omega}_0 = (10^{-3}, 10^{-2}, 10^{-3})$ dává $|\mathbf{L}|^2/(2TI_2) = 1{,}017$; komplementární větev $|\mathbf{L}|^2 < 2 T I_2$ obkružuje osu nejmenšího momentu $\hat{\mathbf{e}}_1$ a plyne záměnou rolí indexů 1 a 3). Dosazením zákonů zachování do Eulerovy rovnice pro $L_y$ a s užitím $dL_y/dt = (I_3 - I_1)(L_z L_x)/(I_3 I_1)$, po algebře:

$$
\boxed{
\left(\frac{dL_y}{dt}\right)^2
\;=\;
\frac{(I_2 - I_1)(I_3 - I_2)}{I_1 I_2^2 I_3}\,
\big[\, \beta^2 - L_y^2 \,\big]\,
\big[\, \gamma^2 - L_y^2 \,\big],
}
$$

kde $\beta^2 = I_2\,(2 T I_3 - |\mathbf{L}|^2)/(I_3 - I_2)$ a $\gamma^2 = I_2\,(|\mathbf{L}|^2 - 2 T I_1)/(I_2 - I_1)$ — nastavené dvěma zachovávajícími se veličinami. Na této větvi $\beta^2 \le \gamma^2$ (rovnost přesně na separatrise $|\mathbf{L}|^2 = 2TI_2$, kde obě rovny $2TI_2$) a $L_y$ osciluje mezi $-\beta$ a $+\beta$: hodnota obratu plyne přímo ze zachování dosazením $L_x = 0$, což dává $L_y^2 = \beta^2$.

Tato ODR je standardním tvarem pro parametrizaci **Jacobiho eliptickými funkcemi**. Uzavřené řešení je (Landau-Lifšic, sv. I §37 Rovn. (37.8)–(37.10); Goldstein 2002 §5.6):

$$
\boxed{
\begin{aligned}
L_x(t) \;&=\; \sqrt{\frac{I_1\,(2TI_3 - |\mathbf{L}|^2)}{I_3 - I_1}}\; \mathrm{cn}\big(\Omega^* t + u_0,\, k\big), \\
L_y(t) \;&=\; \sqrt{\frac{I_2\,(2TI_3 - |\mathbf{L}|^2)}{I_3 - I_2}}\; \mathrm{sn}\big(\Omega^* t + u_0,\, k\big), \\
L_z(t) \;&=\; \sqrt{\frac{I_3\,(|\mathbf{L}|^2 - 2TI_1)}{I_3 - I_1}}\; \mathrm{dn}\big(\Omega^* t + u_0,\, k\big),
\end{aligned}
}
$$

s fází $u_0$ určenou počáteční podmínkou. Kontroly konzistence (přesné, z Jacobiho identit $\mathrm{cn}^2 = 1 - \mathrm{sn}^2$, $\mathrm{dn}^2 = 1 - k^2\mathrm{sn}^2$): $L_x^2 + L_y^2 + L_z^2 = |\mathbf{L}|^2$ a $L_x^2/I_1 + L_y^2/I_2 + L_z^2/I_3 = 2T$ platí identicky pro všechna $t$. Amplituda $L_y$ je $\beta$, ve shodě s argumentem bodu obratu výše; $L_z$ na této větvi nikdy nemizí ($\mathrm{dn} > 0$), což je tvrzení „polhodie obkružuje $\hat{\mathbf{e}}_3$”. Charakteristická frekvence $\Omega^*$ a modul $k$ jsou:

$$
\Omega^*
\;=\;
\sqrt{\frac{(I_3 - I_2)\,(|\mathbf{L}|^2 - 2 T I_1)}{I_1 I_2 I_3}},
\qquad
k^2
\;=\;
\frac{(I_2 - I_1)\,(2 T I_3 - |\mathbf{L}|^2)}{(I_3 - I_2)\,(|\mathbf{L}|^2 - 2 T I_1)}.
$$

(Podrobná algebra je v Landau-Lifšic §37 a Goldstein §5.6 s mírně odlišnými konvencemi; tvar zde sleduje Landau-Lifšic. Implementační poznámka: `scipy.special.ellipk`/`ellipj` přijímají parametr $m = k^2$, nikoli modul $k$.)

> **Oprava (2026-07-11).** Publikovaná podoba tohoto pododdílu (v0.1.0–v0.3.0) nesla chybné zarámované výrazy pro prefaktor redukce, amplitudové koeficienty, $\Omega^*$ a $k^2$ a chybně určovala větev (amplituda $L_y$ byla chybná faktorem ~16 a implikovaná perioda přeskoku o 7–63 % podle konfigurace). Chyba byla nalezena při prvním skutečném provedení křížové kontroly periody přeskoku proti numerické integraci s těsnou tolerancí (forenzní audit 2026-07-11; interní auditní záznam `_audit/forensic-2026-07-11/` ve zdrojovém monorepu): opravené tvary výše souhlasí s numerickým zlatým řešením na $\lesssim 10^{-12}$ v periodě a $8\times10^{-14}$ po prvcích v $(L_x, L_y, L_z)(t)$ napříč konfiguracemi §6.5.5. Regresní pojistka: `tests/test_elliptic_oracle.py`. Jde o čtvrtý incident driftu zarámovaného vzorce (třída KB-043) a první, který dosáhl veřejného zrcadla — právě proto, že křížová kontrola, kterou tento rámeček potřeboval, byla popsána v metadatech fixtur, ale nikdy implementována.

**Perioda eliptického pohybu.** Trajektorie v tělesovém $\mathbf{L}$-prostoru je periodická s periodou:

$$
\boxed{
T_{\mathrm{LL}} \;=\; \frac{4\, K(k)}{\Omega^*},
}
$$

kde $K(k) = \int_0^{\pi/2} d\phi / \sqrt{1 - k^2 \sin^2\phi}$ je úplný eliptický integrál prvního druhu. Toto je **uzavřený výraz pro periodu přeskoku Džanibekova** při konečné amplitudě — analytický obsah OQ-5.5.

**Režim Džanibekova: $|\mathbf{L}|^2$ blízko $2 T I_2$.** V nestabilní rovnováze prostřední osy ($|\mathbf{L}|^2 \to 2 T I_2$, $k \to 1$) perioda diverguje: $K(k) \to \infty$ pro $k \to 1$ ($K(k) \sim \ln(4/\sqrt{1-k^2})$ pro $k \to 1$). Fyzikálně: trajektorie poblíž sedla prostřední osy procházejí *libovolně blízko* něj na časovém měřítku $\sim (1/\sigma) \ln(1/|\delta\Omega|_{\mathrm{initial}})$, ve shodě s linearizovanou rychlostí růstu §6.5.2 (logaritmická závislost na počáteční perturbaci, jak se očekává z dynamiky sedlového bodu). Pro konečné, avšak malé počáteční perturbace $|\delta\Omega_x| \ll |\Omega_2^0|$:

$$
T_{\mathrm{flip}} \;\approx\; \frac{2}{\sigma}\, \ln\!\left(\frac{|\Omega_2^0|}{|\delta\Omega|_{\mathrm{initial}}}\right) \;+\; \mathcal{O}(1),
$$

odpovídajíc linearizované dynamice exponenciálního růstu §6.5.2 v raném režimu. Při velkých počátečních perturbacích ($|\delta\Omega| \sim |\Omega_2^0|$) eliptický modul $k$ klesá hluboko pod 1 a perioda přeskoku asymptoticky směřuje ke konstantě: $T_{\mathrm{flip}} \sim 4 K(k)/\Omega^*$ s $k$ nasyceným na hodnotě určené poměrem setrvačnosti a omezené daleko od 1. **Závislost perioda-přeskoku $\to$ amplituda je tedy $\ln(1/\mathrm{amplituda})$ při malých amplitudách asymptoticky směřující ke konstantě při velkých amplitudách** — falzifikovatelná predikce uzavřeného řešení.

**Geometrická interpretace.** Jacobiho funkce $\mathrm{sn}, \mathrm{cn}, \mathrm{dn}$ parametrizují uzavřenou křivku na průniku polhodie. Pro $k \to 0$ se redukují na $\sin, \cos$ a trajektorie se stává regulární precesí axisymetrického setrvačníku (křížová kontrola: limita $I_1 \to I_2$ dává symetrický setrvačník z §6.3, s $\mathrm{sn}, \mathrm{cn}$ → $\sin, \cos$ a modulem $k \to 0$). Pro $k \to 1$ perioda diverguje a trajektorie tráví většinu času poblíž sedla — režim Džanibekova. Toto jediné uzavřené řešení hladce interpoluje mezi §6.3 volným symetrickým setrvačníkem a §6.5.2 nestabilitou Džanibekova.

### §6.5.4 Podolského ukotvení (CS-jazykové)

**Částečné ukotvení.** Podolský 2018 pokrývá symetrické setrvačníky důkladně (§8.1 volný symetrický setrvačník, §8.2 tlumený, §8.3 Lagrangeův setrvačník), avšak neodvozuje Jacobiho eliptické řešení asymetrického setrvačníku explicitně. Klasickým česko-jazykovým odkazem pro tento materiál je Trkal *Mechanika hmotných bodů a tuhého tělesa* (starší, ale standardní); moderní CS čtenáři mohou též užít Goldstein 2002 v překladu, je-li dostupný. Samotný efekt Džanibekova (pojmenovaný po sovětském kosmonautovi Vladimiru Džanibekovovi, jenž jej pozoroval během své mise na Saljutu-7 roku 1985) je zdokumentován v západní kosmoletové literatuře (Levi 2014 *Classical Mechanics with Calculus of Variations and Optimal Control* §6 — moderní učebnicové pojednání).

### §6.5.5 Numerická tolerance + testovací soubor

| Veličina | Tolerance | Testovací soubor |
|---|---|---|
| Linearizovaná rychlost růstu $\sigma$ vs §6.5.2 analytický vzorec | $< 5\%$ přes lineární režim | `backends/fortran/tests/test_instability_dzhanibekov.f90` (fáze 1 D9 případ 2) + `tests/test_dynamics.py::TestDzhanibekovInstability` (zařazená inženýrská úroveň) |
| Perioda přeskoku o konečné amplitudě vs §6.5.3 zarámovaná $T_{\mathrm{LL}} = 4K(k)/\Omega^*$ | $< 10^{-3}$ relativně; systematicky se zhoršuje poblíž sedlového přiblížení $k \to 1$, jak se $K(k)$ stává citlivým na eliptický modul | (totéž; sken konečné amplitudy fáze 1 D9) |
| Logaritmicko-amplitudové škálování periody přeskoku při malé perturbaci | Shoda sklonu $T_{\mathrm{flip}} \cdot \sigma = 2\ln(|\Omega_2^0|/|\delta\Omega|) + \mathcal{O}(1)$ | (totéž) |
| Zachování $T$, $|\mathbf{L}|^2$ napříč přeskoky | $< 10^{-11}$ relativně jako v testbedu §6.3 | (totéž) |
| Počet dokončených přeskoků přes okno 600 s | Shoda s počtem predikovaným eliptickou funkcí, ±1 | (totéž) |

S uzavřeným řešením z §6.5.3 je celý testbed **třídy A** při lineární i konečné amplitudě. Jacobiho eliptické řešení poskytuje analytická orákula pro *každý* výstup numerické integrace: poloha-na-polhodii, perioda přeskoku, přechod mezi režimy skrze sledování $k$-modulu. Předchozí označení „třída B pouze-numerická“ pro režim o konečné amplitudě je **retirováno** — uzavřené řešení §6.5.3 řeší mezeru analytického orákula.

## §6.6 Testovací případ rozpojení při nekonečné vzdálenosti (třída B)

### §6.6.1 Matematické nastavení + očekávané chování

Pro $|\boldsymbol{\rho}| \to \infty$ člen monopól-monopól $V^{(0)} \to 0$ a oba slapové členy $V_{\mathrm{tidal},A}, V_{\mathrm{tidal},B} \to 0$ jako $1/|\boldsymbol{\rho}|^3$. Celkový vzájemný potenciál mizí rychleji než jakékoli inverzně-mocninné fyzikální měřítko. Každé těleso se stává izolovaným bezsilovým setrvačníkem: $\boldsymbol{\tau}_A, \boldsymbol{\tau}_B \to \mathbf{0}$ (úměrné $1/|\boldsymbol{\rho}|^3$), $V \to 0$ a Lagrangián se faktorizuje jako $L = L_A + L_B$, kde každé $L_X = T_X$ je kineticky-pouze Lagrangián izolovaného tuhého tělesa. Soustava se zcela rozpojí.

V praxi: umísti tělesa do velké vzdálenosti $|\boldsymbol{\rho}| = R_{\mathrm{far}}$, zvol počáteční tělesové úhlové rychlosti tak, aby rotační dynamika každého tělesa byla nezávislou bezsilovou Eulerovou precesí (kap. 06 §6.3) příslušné třídy symetrie, a ověř numericky, že křížová vazba (zbytkový moment síly z gravitačního gradientu + zpětná síla při této velké, ale konečné $R_{\mathrm{far}}$) zůstává pod prahem předepsaným asymptotickým škálováním.

### §6.6.2 Očekávané škálování + tolerance

Při vzdálenosti $R_{\mathrm{far}}$ je velikost zbytkového momentu síly z gravitačního gradientu $|\boldsymbol{\tau}_A| \sim 3Gm_B \cdot |\mathbf{I}_A - I_0 \mathbb{I}|/R_{\mathrm{far}}^3$. Volbou $R_{\mathrm{far}} = 10^3 \cdot \ell_A$, kde $\ell_A$ je rozměr tělesa, je bezrozměrná zbytková vazba $\sim 10^{-9}$ vůči bezsilovému spinovému momentu hybnosti. Přes jednu orbitální periodu $T = R_{\mathrm{far}}^{3/2}/(GM)^{1/2}$ by akumulovaný efekt na individuální moment hybnosti každého tělesa měl zůstat pod $10^{-12}$.

### §6.6.3 Testovací soubor

| Veličina | Tolerance | Testovací soubor |
|---|---|---|
| Křížová vazba mezi momenty hybnosti tělesa-$A$ a tělesa-$B$ | $< 10^{-12}$ relativně na orbitální periodu | `tests/test_dynamics.py::TestInfiniteSeparationDecoupling` (plánováno, dle `_specs/...` §3.2) |
| Stabilita individuálního $|\mathbf{L}_X|$ každého tělesa | Tolerance třídy A §6.3 | (totéž) |
| Stabilita kinetické energie $T_X$ každého tělesa | Tolerance třídy A §6.3 | (totéž) |

Tento testbed validuje **poctivost multipolového useknutí** kap. 03 §3.4.3 zarámovaného: členy vyššího multipolového řádu by měly klesat správně s $1/r^n$ a žádná nefyzikální vazba při velkém $r$ by neměla vzniknout z numerických artefaktů.

## §6.7 Bitová identita mezi implementacemi (třída B procedurální — doktrína §4.4 C2)

### §6.7.1 Disciplína

Když existují dvě nezávislé implementace zarámované sady rovnic kap. 04 §4.7 — nyní: (a) `dynamics.py` Python první verze na inženýrské úrovni, (b) `backends/fortran/src/jwst_l2_dynamics.f90` fortranový lege-artis backend ve fázi 1 — výstupy mezi implementacemi musí souhlasit při doktrínou předepsané toleranci. Podle kritéria C2 z `_config/LEGE-ARTIS-LANGUAGE-DOCTRINE-v0.1.md` §4.4 (aplikovaného v doméně analytické mechaniky) je tolerancí **ULP × √N** dle Higham §4.2 Option G, kde ULP je jednotka na posledním místě při velikosti srovnávané veličiny a $N$ je počet kumulativních aritmetických operací (kroků RK4 v našem případě).

Pro konfiguraci první verze: $N \approx 12\,000$ kroků RK4 přes 600 s při $dt = 0{,}05$ s; ULP při $|\boldsymbol{\rho}| \sim 10$ m je $\sim 10^{-15}\, |\boldsymbol{\rho}|$; ULP × √N $\sim 10^{-13}$ relativně. Toto je **brána** pro srovnání mezi implementacemi.

### §6.7.2 Implementace testbedu

Srovnání mezi implementacemi běží napříč těmito výstupy fáze 1 (dle `_handoffs/fortran-phase-1/SONNET-BRIEF.md`):

- **D4 — křížový test momentu síly z gravitačního gradientu.** Jediná konfigurace; srovnej fortranem-vypočtené $\boldsymbol{\tau}_X^{\mathrm{body}}$ vs pythonem-vypočtené při identických vstupech. Tvrdá brána při relativní toleranci ULP × √N.
- **D8 — křížový test trajektorie symetrického setrvačníku.** Jediná 12\,000kroková RK4 trajektorie konfigurace testbedu §6.3; srovnej plnou časovou řadu $(t, \boldsymbol{\Omega}_X(t))$. Tvrdá brána při relativní toleranci ULP × √N.

Po uzavření fáze 1 (tag `jwst-l2-fortran-v0.1.0-phase-1-green`) se disciplína mezi implementacemi rozšiřuje na plnou sadu rovnic, s fází 2 pokrývající testy librace + Džanibekova.

### §6.7.3 Co testuje

Toto je analogem *knihovnou-podloženého vs ručně-psaného FFT* z projektu fourier (kap. 03 §3.6 + záznamy projektového vlastníka JWST-L2 §3.3) v doméně analytické mechaniky: když 2+ implementace souhlasí bit po bitu, pravděpodobnost, že obě sdílejí tutéž algebraicko-implementační chybu, dramaticky klesá. Bitová identita mezi implementacemi nechytá chyby v samotných rovnicích kanonické úrovně (obě implementace by mohly implementovat tutéž chybnou rovnici); pro to jsou relevantní pojistkou testbedy třídy A vůči analytickým řešením kap. 06 §6.2–§6.5. C2 chytá **implementační úlet** mezi kanonickými rovnicemi a jejich numerickou realizací.

## §6.8 Co je stanoveno na konci této kapitoly

Čtenář, jenž prošel §6.1–§6.7, má:

- **Tvrzení disciplíny** (§6.1): uzavřená analytická řešení třídy A jako orákula pro integrátor; konfigurace okrajových případů třídy B jako validátory domény modelu. Obojí nutné pro poctivou validaci.
- **Základní testbedy třídy A:**
  - §6.2 Keplerova limita (analytická Keplerova rovnice v limitě $\mathbf{I}_X \to \mathbf{0}$; tolerance $10^{-6}$ relativně na orbitální periodu).
  - §6.3 Volný symetrický setrvačník (Podolský §8.1 Rovn. 8.1–8.2; obnovení precesní frekvence v rámci 0,2 FFT binů; implementováno v první verzi s $|dE/E_0| = 6{,}2 \times 10^{-12}$).
  - §6.4 Librace gravitačního gradientu linearizovaná (Hughes §8.4 + Podolský §8.3 metoda efektivního potenciálu; librační perioda v rámci $10^{-4}$ přes 10 orbit; testbed fáze 2).
  - §6.5 Nestabilita Džanibekova linearizovaná (analytická rychlost růstu $\sigma$; fáze 1 D9 případ 2).
- **Základní testbedy třídy B:**
  - §6.5 Džanibekov plně-nelineární (dynamika přeskoků, perioda závislá na amplitudě).
  - §6.6 Rozpojení při nekonečné vzdálenosti (křížová vazba $< 10^{-12}$ na orbitální periodu).
  - §6.7 Bitová identita mezi implementacemi (doktrína §4.4 C2 při ULP × √N).
- **Mapu ukazatelů** z každé specializace třídy symetrie v kap. 05 na její odpovídající testbed v této kapitole.

**Co dosud stanoveno NENÍ ve v0.1:**

- **Testbedy omezené soustavy tří těles v soustavě L2** (linearizované vlastní hodnoty stability u L2; uzavřené přiblížení halo orbity; Ljapunovovy orbity). Rozsah kanonické úrovně v0.2.
- **Symplektický integrátor na Lieově grupě** jako alternativní implementace k validaci C2 brány mezi implementacemi. Současný C2 testbed srovnává Python ↔ Fortran (dvě RK4 implementace týchž rovnic); port symplektického integrátoru by umožnil bohatší C2 křížovou validaci. Rozsah inženýrské úrovně v0.2 (kap. 04 §4.9 OQ-4.6).
- **Plné Ljapunovovo spektrum** asymetrického tělesa v poli gravitačního gradientu (relevantní zobecnění testbedu §6.5 Džanibekova na případ perturbovaný gravitačním gradientem, s kvantitativní predikcí Ljapunovova exponentu). Rozsah kanonické úrovně v0.2 (kap. 04 §4.9 OQ-4.8).
- **Testbed relativistických 1PN-korekcí** pro režim $v/c \to \mathrm{nemalé}$. Explicitně mimo rozsah v0.1 (kap. 00 §6); mohl by se stát záznamem třídy A ve v0.3+, pokud se rozsah kdy rozšíří.

---

## §6.9 Otevřené otázky v této kapitole

### Dopředu orientované (rozsah v0.2)

- **OQ-6.1.** Implementovat a validovat testbedy **omezené soustavy tří těles u L2**: linearizované vlastní hodnoty stability u L2 (Szebehely 1967 §4; Murray & Dermott 1999 §3); uzavřené přiblížení halo orbity (Richardson 1980 rozvoj třetího řádu; Howell 1984); numerické pokračování rodiny Ljapunovových orbit. Každý je testbedem třídy A nebo B pro kanonickou úroveň v0.2 (jež přidává rotující soustavu Země-Slunce + centrifugální + Coriolisovy pseudo-síly).
- **OQ-6.2.** Implementovat **posouzení dlouhohorizontového úletu RK4** (`_specs/...` §3.7): integrovat testbed §6.3 přes $T = 10^4 T_{\mathrm{orbit}}$ a kvantifikovat sekulární úlet energie, potvrzujíc škálování $O(dt^4 \cdot T)$ pro nesymplektický RK4. Toto je empirickou základnou, vůči níž je symplektický integrátor na Lieově grupě (kap. 04 §4.9 OQ-4.6) srovnáván.
- **OQ-6.3.** Implementovat **audit řádu přesnosti konvergence RK4** (`_specs/...` §3.6): integrovat testbed §6.3 při $dt \in \{0{,}5, 0{,}1, 0{,}05, 0{,}01, 0{,}005, 0{,}001\}$ s a potvrdit, že log-log proložení $|dE|$ vs $dt$ má sklon $4{,}0 \pm 0{,}1$ v konvergentním režimu. Analog doktríny §4.4 C4 pro integraci ODR.
- **OQ-6.4** ***Vyřešeno*** — v §6.4.3 (tato kapitola): uzavřené řešení Lagrangeova setrvačníku + librace gravitačního gradientu metodou efektivního potenciálu $V_{\mathrm{ef}}(\vartheta)$, s prvními integrály cyklických souřadnic $L_\psi, L_\varphi$, zarámovaným efektivním potenciálem, uzavřenou kvadraturou pro $\vartheta(t)$ a třemi precesními režimy (a)/(b)/(c). Falzifikovatelné škálování amplituda-perioda $T_{\mathrm{lib}}(\vartheta_{\mathrm{amp}}) = (2/\pi)K(k)\cdot T_{\mathrm{lib}}^{\mathrm{linear}}$. Páruje se s řešením kap. 05 §5.8 OQ-5.6.
- **OQ-6.5** ***Vyřešeno*** — v §6.5.3 (tato kapitola): uzavřená parametrizace Jacobiho eliptickými funkcemi $(L_x, L_y, L_z) = (\alpha\sqrt{\cdot}\,\mathrm{cn}, \alpha\,\mathrm{sn}, \alpha\sqrt{\cdot}\,\mathrm{dn})$ s eliptickým modulem $k$ a charakteristickou frekvencí $\Omega^*$ určenými $T$ a $|\mathbf{L}|^2$. Zarámovaná perioda $T_{\mathrm{LL}} = 4K(k)/\Omega^*$ → perioda přeskoku Džanibekova jako falzifikovatelná na amplitudě závislá veličina. Páruje se s řešením kap. 05 §5.8 OQ-5.5.

### Procedurální / inženýrská úroveň

- **OQ-6.6.** Implementovat C2 bránu bitové identity mezi implementacemi (§6.7) jako CI workflow na běžci GitHub Actions pro oba Python a Fortran backendy. Vzorováno podle CI bitové-identity-mezi-backendy v0.2.x projektu `lege-artis/fourier` (dle fourier doktríny §4.4 C2 briefu rozšíření testbedu).

---

> **Doplňující anglojazyčné odkazy.** Pro testbedy analytické mechaniky: Landau-Lifšic, sv. I §13–§15 (Kepler), §35 (symetrický setrvačník), §37 (asymetrický setrvačník); Goldstein 2002 §5.4–§5.7 (mechanika tuhého tělesa včetně Džanibekova + Lagrangeova setrvačníku); Hughes 2004 §3–§5 + §8 (dynamika postoje kosmické lodě + stabilizace gravitačním gradientem). Pro moderní perspektivu Lieových grup / variační: Marsden & Ratiu 1999 §13 (Euler-Poincaré); Holm-Schmah-Stoica 2009 §6 (tuhé těleso $SO(3)$ jako pracovní příklad). Pro rámec symplektického / Lie-grupového integrátoru (párování kap. 04 §4.9 OQ-4.6): Hairer-Lubich-Wanner 2006 §VII.5 (symplektický Runge-Kutta na Lieových grupách). Pro rozšíření omezené soustavy tří těles u L2 (v0.2): Szebehely 1967; Murray & Dermott 1999 §3.
>
> **Hlavní odkaz.** Podolský 2018 Kapitola 8 *Aplikace: setrvačníky* poskytuje přímá CS-jazyková orákula analytických limitů pro tři z šesti základních testbedů: §8.1 (Rovn. 8.1–8.2) pro §6.3 volný symetrický setrvačník; §8.2 tlumená varianta; §8.3 (Rovn. 8.10–8.15) pro §6.4 libraci gravitačního gradientu přes Lagrangeův setrvačník. Testbed §6.2 Kepler má klasické anchory v CS literatuře mimo Podolského (překlad Landau-Lifšic; starší Trkal). Asymetrický setrvačník §6.5 Džanibekov *není* přímo pokryt v Podolském — moderní CS čtenáři by měli párovat s Goldstein 2002 §5.6.

---

**Následující:** v0.2 kanonické úrovně (rozšíření omezené soustavy tří těles u L2; obsah kapitoly o symplektickém integrátoru na Lieově grupě; plná nelineární Lagrangeova-setrvačníková + Jacobiho-eliptická asymetricko-setrvačníková analytická řešení).

**Konec 06-analyticke-limity-a-validacni-brany.md (v0.1, 3. kolo).**
