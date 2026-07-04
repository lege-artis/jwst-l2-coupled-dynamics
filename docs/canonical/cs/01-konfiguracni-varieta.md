# Kapitola 1 — Konfigurační varieta a variační princip

> **Předpokládané čtení.** Kapitola 00 (úvod). Znalost elementárního vícerozměrného kalkulu a pojmu hladké variety. Žádná předchozí zkušenost se symplektickou geometrií není nutná; struktury jsou zaváděny podle potřeby.
>
> **Stav.** v0.1 — prvních 3–5 stran z ~15stránkové zamýšlené kapitoly. Následující oddíly budou doplněny v samostatných autorských posezeních Opus.

> **Poznámka ke způsobu a předpokladům odvození.** Tato CS kanonická kapitola (v0.1) je odvozena z anglické kanonické kapitoly `docs/canonical/en/01-configuration-manifold.md` (v0.1). Terminologie je ukotvena v Podolský 2018 *Teoretická mechanika ve třech knihách* (Matfyzpress). EN verze zůstává autoritativním pramenem při budoucích aktualizacích.

---

## §1.1 Co tato kapitola dělá

Tato kapitola připravuje geometrickou scénu pro vše, co následuje. Identifikujeme **konfigurační prostor** soustavy dvou tuhých těles — množinu všech možných *poloh a orientací*, jež dvojice těles může zaujmout — a vybavíme jej strukturou potřebnou k zapsání Lagrangeovy mechaniky na něm: tečným svazkem, hladkou Lagrangeovou funkcí a variačním principem. Na konci kapitoly budou Eulerovy–Lagrangeovy pohybové rovnice naší soustavy připravené k odvození (kapitola 04) a otázka „v jaké konkrétní soustavě máme počítat?“ se stane vědomou volbou, nikoli předpokladem.

Pojetí je **nezávislé na volbě soustavy**: popisujeme soustavu jazykem, jenž neprivileguje žádnou konkrétní vztažnou soustavu. Když nakonec soustavu zvolíme (od kapitoly 05), činíme tak proto, že rovnice jsou tam nejjednodušší, nikoli proto, že jsme zapomněli, že máme volbu.

## §1.2 Konfigurace jednoho tuhého tělesa

Než popíšeme dvě provázaná tělesa, popíšeme jedno.

**Tuhé těleso** je soubor hmotných bodů $\{m_\alpha\}$ v polohách $\{\mathbf{r}_\alpha\}$ takový, že všechny párové vzdálenosti $|\mathbf{r}_\alpha - \mathbf{r}_\beta|$ jsou v čase konstantní. Podle klasické věty (Goldstein 2002 §4.1; Arnold 1989 §3.27) je konfigurace tuhého tělesa v trojrozměrném euklidovském prostoru zcela určena:

1. Polohou jednoho zvoleného referenčního bodu (běžně těžiště) v jisté inerciální soustavě: bodem $\mathbf{R} \in \mathbb{R}^3$.
2. Orientací tělesa vůči této inerciální soustavě: rotací $q \in SO(3)$.

Zde $SO(3)$ značí **speciální ortogonální grupu** ve třech dimenzích — množinu reálných matic $3 \times 3$ $\mathbf{R}_q$ splňujících $\mathbf{R}_q^T \mathbf{R}_q = \mathbf{I}$ a $\det(\mathbf{R}_q) = +1$. Geometricky je $SO(3)$ množinou rotací trojprostoru, jež zachovávají orientaci (pravotočivost). Jako hladká varieta má $SO(3)$ dimenzi tři; jednou běžnou parametrizací jsou Eulerovy úhly $(\varphi, \vartheta, \psi)$, jinou jednotkové kvaterniony $q \in S^3 \subset \mathbb{R}^4$ s identifikací $q \sim -q$ (dvojnásobné nakrytí; viz Marsden & Ratiu 1999 §9.1). Kód inženýrské úrovně z první verze používá jednotkové kvaterniony; odvození kanonické úrovně je na parametrizaci nezávislé.

Konfigurace jednoho tuhého tělesa je tedy prvkem:

$$
SE(3) = \mathbb{R}^3 \times SO(3),
$$

**speciální euklidovské grupy** tuhých pohybů v trojprostoru. Jde o šestirozměrnou Lieovu grupu. První faktor popisuje translaci, druhý rotaci.

> **Poznámka.** Struktura polopřímého součinu $SE(3) = \mathbb{R}^3 \rtimes SO(3)$ je nezbytná, když probíráme, jak se rotace a translace skládají (rotovat posunutý objekt není totéž jako posunout rotovaný). Tuto strukturu budeme potřebovat explicitně v kapitole 04 při probírání redukce podle symetrií. Prozatím postačí $SE(3) \cong \mathbb{R}^3 \times SO(3)$ jako varieta.

## §1.3 Konfigurace dvou tuhých těles

Pro dvě tuhá tělesa $A$ a $B$, každé nezávisle umístěné a orientované v trojprostoru, je konfigurační prostor kartézským součinem:

$$
Q \;=\; SE(3)_A \times SE(3)_B \;=\; \big(\mathbb{R}^3 \times SO(3)\big)_A \;\times\; \big(\mathbb{R}^3 \times SO(3)\big)_B.
$$

Bod v $Q$ je čtveřice:

$$
q \;=\; (\mathbf{R}_A,\ q_A,\ \mathbf{R}_B,\ q_B), \qquad \mathbf{R}_A, \mathbf{R}_B \in \mathbb{R}^3, \quad q_A, q_B \in SO(3).
$$

Dimenze $Q$ je $6 + 6 = 12$. To je očekávaný počet: tři polohy pro těžiště tělesa $A$, tři úhly pro orientaci tělesa $A$, totéž pro $B$. Kód inženýrské úrovně z první verze (`dynamics.py`, třída `BodyState`) nese 13 složek na těleso (3 poloha + 3 rychlost + 4 kvaternion + 3 úhlová rychlost = 13), což je počet *fázového prostoru*, dvojnásobek počtu konfiguračního plus přepočet kvaternion-versus-rotační-matice.

### §1.3.1 Bezprostřední symetrie — a co s ní uděláme

Konfigurační prostor $Q$ má na sobě přirozenou akci samotné grupy $SE(3)$: je-li $g = (\mathbf{a}, h) \in SE(3)$ globální tuhý pohyb (translace o $\mathbf{a}$, rotace o $h$), pak

$$
g \cdot q \;=\; \big(h \mathbf{R}_A + \mathbf{a},\ h \cdot q_A,\ h \mathbf{R}_B + \mathbf{a},\ h \cdot q_B\big)
$$

popisuje tutéž fyzikální konfiguraci jako $q$, jen pohlíženou z posunuté a otočené vztažné soustavy.

Toto je formálním vyjádřením tvrzení, že **fyzika nezávisí na volbě inerciální vztažné soustavy**. Akce grupy $SE(3)$ na $Q$ převádí konfiguraci na jinou konfiguraci, jež je fyzikálně ekvivalentní.

Faktorprostor $Q / SE(3)$ je prostorem fyzikálně odlišných konfigurací. Jde o šestirozměrnou varietu (dvanáct minus šest dimenzí grupy symetrie). Nebudeme pracovat přímo ve faktorprostoru — pro výpočetní účely je jednodušší zvolit soustavu, jež symetrii fixuje — avšak existence tohoto faktorprostoru je tím, co nám později (kapitola 02) umožní zvolit **těžišťovou soustavu** bez ztráty obecnosti. Těžišťová soustava fixuje tři z šesti translačně/rotačních parametrů symetrie; zbylé tři jsou fixovány srovnáním os s vektorem relativní polohy nebo volbou orientace pevné vůči tělesu. Žádná z těchto voleb nemění fyziku.

> **Konvence přijatá od tohoto bodu.** Pokud není uvedeno jinak, pracujeme v **inerciální soustavě** s počátkem v těžišti soustavy dvou těles. Translační invariance pohybových rovnic implikuje, že toto těžiště je samo inerciální. Celková lineární hybnost $\mathbf{P}_{\text{total}} = m_A \dot{\mathbf{R}}_A + m_B \dot{\mathbf{R}}_B$ je konstantní; volbou soustavy ji klademe rovnu nule. Šest zachovávajících se veličin (poloha $\mathbf{R}_{\text{CoM}}$ a hodnota $\mathbf{P}_{\text{total}}$) je tím trivializováno. Zbylá dynamika žije na šestirozměrném redukovaném konfiguračním prostoru, výhodně parametrizovaném jako $(\Delta \mathbf{R},\ q_A,\ q_B)$, kde $\Delta \mathbf{R} = \mathbf{R}_B - \mathbf{R}_A$ je vektor relativní polohy a $q_A$, $q_B$ zůstávají jako dříve. Tuto redukci budeme volně užívat od kapitoly 02.

## §1.4 Rychlosti — tečný svazek

Lagrangeova mechanika vyžaduje nejen polohy, ale i rychlosti. Množina všech dvojic (poloha, rychlost) v dané konfiguraci je **tečný prostor** ke konfigurační varietě v této konfiguraci; disjunktní sjednocení tečných prostorů přes všechny konfigurace je **tečný svazek** $TQ$.

Pro naši soustavu má tečný svazek dimenzi $2 \cdot 12 = 24$. Bod v $TQ$ je konfigurace $q$ spolu s vektorem rychlosti:

$$
(q,\ \dot{q}) \;=\; \big( (\mathbf{R}_A, q_A, \mathbf{R}_B, q_B),\ (\dot{\mathbf{R}}_A, \dot{q}_A, \dot{\mathbf{R}}_B, \dot{q}_B) \big).
$$

Složky translační rychlosti $\dot{\mathbf{R}}_A, \dot{\mathbf{R}}_B \in \mathbb{R}^3$ jsou bezproblémové — žijí v tečném prostoru $\mathbb{R}^3$, jímž je $\mathbb{R}^3$ sám. Složky rotační rychlosti $\dot{q}_A, \dot{q}_B$ vyžadují opatrnost, neboť $SO(3)$ není vektorový prostor; jeho tečný prostor v bodě *není* totéž co varieta sama.

Tečný prostor $T_{q_A} SO(3)$ v orientaci $q_A \in SO(3)$ je trojrozměrný vektorový prostor. Lze jej ztotožnit s Lieovou algebrou $\mathfrak{so}(3)$ antisymetrických matic $3 \times 3$, jež je zase přirozeně izomorfní s $\mathbb{R}^3$ skrze **operátor stříšky (hat)**:

$$
\hat{}\colon \mathbb{R}^3 \to \mathfrak{so}(3), \qquad
\boldsymbol{\omega} = (\omega_1, \omega_2, \omega_3) \;\mapsto\;
\hat{\boldsymbol{\omega}} \;=\;
\begin{pmatrix}
0 & -\omega_3 & \omega_2 \\
\omega_3 & 0 & -\omega_1 \\
-\omega_2 & \omega_1 & 0
\end{pmatrix}.
$$

Toto zobrazení je izomorfismem Lieových algeber (vektorový součin na $\mathbb{R}^3$ se stává komutátorem matic na $\mathfrak{so}(3)$): $\widehat{\mathbf{a} \times \mathbf{b}} = [\hat{\mathbf{a}}, \hat{\mathbf{b}}]$.

Fyzikální obsah operátoru stříšky: je-li dána rotační matice $\mathbf{R}_q(t)$ závislá na čase, její časová derivace je $\dot{\mathbf{R}}_q = \hat{\boldsymbol{\omega}}\, \mathbf{R}_q$, kde $\boldsymbol{\omega} \in \mathbb{R}^3$ je **úhlová rychlost v inerciální soustavě**. Ekvivalentně $\dot{\mathbf{R}}_q = \mathbf{R}_q\, \hat{\boldsymbol{\Omega}}$, kde $\boldsymbol{\Omega} = \mathbf{R}_q^T \boldsymbol{\omega}$ je **úhlová rychlost v tělesové soustavě**. Obě reprezentace popisují tutéž fyzikální rychlost rotace; volba mezi nimi je otázkou, jaká soustava je nejvýhodnější.

> **Křížový odkaz na kód.** `dynamics.py::q_kinematics_matrix` z první verze počítá $\dot{q}$ pro reprezentaci $\mathbf{R}_q$ jednotkovým kvaternionem; vzorec je $\dot{q} = \tfrac{1}{2}\, q \cdot (0, \boldsymbol{\Omega})$ v kvaternionové algebře. Jde o tělesovou ($\boldsymbol{\Omega}$) verzi výše uvedeného vztahu. Odvození: vezmi časovou derivaci $\mathbf{R}_q^T \mathbf{R}_q = \mathbf{I}$, čímž získáš $\dot{\mathbf{R}}_q^T \mathbf{R}_q + \mathbf{R}_q^T \dot{\mathbf{R}}_q = 0$, tudíž $\mathbf{R}_q^T \dot{\mathbf{R}}_q$ je antisymetrická a rovná se $\hat{\boldsymbol{\Omega}}$ pro jisté $\boldsymbol{\Omega} \in \mathbb{R}^3$. Dosaď kvaternionovou parametrizaci a zjednoduš.

Prakticky to znamená, že **rychlosti na naší 12rozměrné konfigurační varietě jsou popsány 12 čísly**: třemi složkami translační rychlosti pro těžiště každého tělesa a třemi složkami úhlové rychlosti pro rotaci každého tělesa (buď v inerciální, nebo v tělesové soustavě, dle volby). Tečný svazek $TQ$ je dvanáctirozměrný v konfiguraci plus dvanáctirozměrný v rychlosti, celkem dvacet čtyři.

## §1.5 Lagrangián a variační princip

Máme konfigurační varietu $Q$ a její tečný svazek $TQ$. **Lagrangeova mechanika** sestává z:

1. Hladké funkce $L \colon TQ \times \mathbb{R} \to \mathbb{R}$ — **Lagrangiánu** — přiřazující reálné číslo každé trojici (konfigurace, rychlost, čas). Pro naši autonomní (na čase nezávislou) soustavu nemá $L = L(q, \dot{q})$ explicitní časovou závislost.
2. Variačního principu: skutečná fyzikální trajektorie mezi dvěma konfiguracemi $q(t_1)$ a $q(t_2)$ je ta, jež činí **akční funkcionál**

$$
S[q] \;=\; \int_{t_1}^{t_2} L(q(t), \dot{q}(t))\ dt
$$

extremálním vzhledem k malým variacím dráhy $q(t)$, jež se v koncových bodech anulují. Toto je **Hamiltonův princip** (Goldstein §2.1; Thorne & Blandford sv. I §7.2; Feynman sv. II kap. 19; Arnold §3.13). Pro mechanické soustavy s kinetickou energií $T$ a potenciální energií $V$ má Lagrangián tvar $L = T - V$.

> **Proč variační princip?** Tři odpovědi, ve vzestupném pořadí hloubky:
> 1. **Operačně**: variační formulace činí symetrie a zákony zachování transparentními skrze Noetherové větu (§1.6 níže).
> 2. **Geometricky**: akční funkcionál a jeho extrémy jsou objekty nezávislé na souřadnicích; výsledné pohybové rovnice mají zaručeně tentýž tvar v jakékoli volbě zobecněných souřadnic.
> 3. **Foundačně**: rozšíření klasické mechaniky na kvantovou je přirozenější v lagrangeovsko-akční formulaci (Feynmanovy dráhové integrály) než v newtonovské formulaci síla-a-zrychlení. Tuto perspektivu pro zdejší kanonickou úroveň nepotřebujeme, ale její existence je silnou kontrolou zdravého rozumu, že jde o správný výchozí bod.

Podmínka extremality $\delta S = 0$ pro variace $\delta q(t)$ anulující se v $t_1, t_2$ dává **Eulerovy–Lagrangeovy rovnice**:

$$
\frac{d}{dt} \frac{\partial L}{\partial \dot{q}^i} \;-\; \frac{\partial L}{\partial q^i} \;=\; 0, \qquad i = 1, \ldots, \dim Q.
$$

Zde $q^i$ jsou lokální souřadnice na $Q$ a $\dot{q}^i$ jejich časové derivace. Odvození této rovnice z $\delta S = 0$ je standardní (Goldstein §2.3; Arnold §3.14); vynecháváme jej, neboť do něj nevstupuje žádná struktura specifická pro naši soustavu dvou těles.

Pro naši konkrétní soustavu provádějí konstrukci podrobně kapitoly 02 (kinetická energie na $SE(3) \times SE(3)$), 03 (vzájemný gravitační potenciál) a 04 (sestavení $L = T - V$ a aplikace Eulerova–Lagrangeova operátoru).

## §1.6 Noetherové věta — symetrie se stávají zákony zachování

Spojitá symetrie Lagrangiánu — tj. jednoparametrická rodina transformací $q \mapsto q_\epsilon$ taková, že $L(q_\epsilon, \dot{q}_\epsilon) = L(q, \dot{q})$ pro všechna $\epsilon$ blízká nule — implikuje zachovávající se veličinu podél fyzikálních trajektorií. Toto je **Noetherové věta** (Goldstein §13.7; Arnold §A.7; Marsden & Ratiu §11.2).

Pro naši soustavu dvou těles jsou bezprostředními symetriemi Lagrangiánu (jež formálně ověříme v kapitole 04):

| Symetrie | Zachovávající se veličina | Diagnostika zachování v inženýrské úrovni |
|---|---|---|
| Časová translace $t \to t + \epsilon$ | Celková energie $E = T + V$ | `dynamics.py::total_kinetic_energy + total_potential_energy`; první verze dosahuje $\|dE/E_0\| = 6{,}2 \times 10^{-12}$ |
| Prostorová translace $\mathbf{R}_A, \mathbf{R}_B \to \mathbf{R}_A + \boldsymbol{a}, \mathbf{R}_B + \boldsymbol{a}$ | Celková lineární hybnost $\mathbf{P}$ | První verze: $\mathbf{P}$ zůstává nulová na strojovou přesnost, neboť těžišťová soustava je vynucena při inicializaci |
| Prostorová rotace o $h \in SO(3)$ působící na obě tělesa | Celkový moment hybnosti $\mathbf{L}$ | `dynamics.py::total_angular_momentum`; první verze dosahuje $\|d\mathbf{L}/\mathbf{L}_0\| = 2{,}1 \times 10^{-10}$ |

Třetí řádek je strukturálně důležitý: celkový moment hybnosti se zachovává právě proto, že Lagrangián je invariantní vůči současné globální rotaci obou těles. Toto *není* totéž jako zachování momentu hybnosti *každého jednotlivého tělesa* (moment hybnosti každého tělesa se vyměňuje s druhým skrze moment síly z gravitačního gradientu). Symetrie, jež dává zachování celkového $\mathbf{L}$, je rotační invariance *vzájemné* potenciální energie.

> **Mezi­úrovňová poznámka.** Noetherové věta je tím, co *slibuje*, že diagnostiky zachování z první verze budou fungovat na strojovou přesnost v nepřítomnosti numerické chyby. Variační struktura podkladových pohybových rovnic to zaručuje. RK4 (integrátor inženýrské úrovně) je nesymplektický, a vnáší tedy malou chybu zachování energie a momentu hybnosti na krok; pozorované škály $10^{-12}$ / $10^{-10}$ jsou kumulativní integrační chybou za 600 sekund. Symplektický integrátor na Lieově grupě (Hairer–Lubich–Wanner 2006 §VII.5) by posunul zachování energie blíže k zaokrouhlovací chybě; toto probíráme v kapitole 06.

## §1.7 Co je stanoveno na konci této kapitoly

Čtenář, jenž prošel §1.1–§1.6, má:

- Konfigurační prostor $Q = SE(3) \times SE(3)$ naší soustavy, s explicitní dimenzí 12.
- Tečný svazek $TQ$ dimenze 24, s rychlostmi parametrizovanými translačními rychlostmi $\dot{\mathbf{R}}_A, \dot{\mathbf{R}}_B$ a úhlovými rychlostmi $\boldsymbol{\omega}_A, \boldsymbol{\omega}_B$ (buď v inerciální, nebo v tělesové soustavě, dle naší volby).
- Variační princip vybírající fyzikální trajektorii mezi danými koncovými body: $\delta S = 0$ pro $S = \int L\, dt$.
- Eulerovy–Lagrangeovy rovnice jako lokální tvar variačního principu.
- Noetherové větu spojující symetrie $L$ se zákony zachování, připravenou k aplikaci, jakmile je $L$ specifikováno.

**Co dosud stanoveno NENÍ:**

- Kinetická energie $T$ explicitně jako funkce na $TQ$. To vyžaduje tenzory setrvačnosti obou těles. **Kapitola 02.**
- Potenciální energie $V$ explicitně. To vyžaduje vzájemný gravitační potenciál jako funkci konfigurace. **Kapitola 03.**
- Eulerovy–Lagrangeovy rovnice v explicitním souřadnicovém tvaru pro náš konkrétní Lagrangián. **Kapitola 04.**
- Eulerovy rovnice v tělesové soustavě jako specializace Euler-Poincarého redukcí. **Kapitola 05.**
- Validace vůči analytickým limitům. **Kapitola 06.**

> **Hlavní odkaz.** Pro elementární souřadnicový aparát (ortogonální transformační matice + Eulerovy úhly) je hlavním CS-jazykovým ukotvením Podolský 2018 Kapitola 6 *Kinematika tuhého tělesa*: §6.1 *Vektory a tenzory* (tenzorová algebra v euklidovském prostoru, klasické vs moderní definice), §6.2 *Relativita otáčivého pohybu* (báze inerciální vs pevná v tělese $\mathbf{e}'_i = A_{ik}(t)\mathbf{e}_k$ s ortogonální $A$), §6.4 *Eulerovy úhly* (intrinzická-ZXZ parametrizace $A = BCD$ s Rovn. (6.14)–(6.16) pro tři elementární rotační matice). Podolský *neabstrahuje* na úroveň Lieových grup — zachází s rotacemi jako s ortogonálními maticemi $3 \times 3$ parametrizovanými Eulerovými úhly nebo jednotkovými kvaterniony, nikoli jako s prvky hladké variety $SO(3)$ s její strukturou tečného svazku a levo-invariantních vektorových polí.
>
> **Doplňující anglojazyčné odkazy.** Pro abstrakci na úroveň Lieových grup (užitou v této kapitole a využitou v kapitole 02 §2.7.1 a kapitole 04 §4.4), jež je moderní Marsden-Ratiu materiál nepřítomný v Podolském: Arnold 1989 *Mathematical Methods of Classical Mechanics* §§1.4, 4.A (konfigurační variety + lagrangeovské soustavy na $TM$); Marsden & Ratiu 1999 *Introduction to Mechanics and Symmetry* §§1, 7, 9.1 ($SO(3)$ a $SE(3)$ jako Lieovy grupy se strukturou tečného svazku); Sussman & Wisdom 2015 *Structure and Interpretation of Classical Mechanics* §1 (výpočetní vykreslení téhož materiálu); pro $SO(3)$ jako Lieovu grupu specificky Holm-Schmah-Stoica 2009 §6 (příprava na Euler-Poincarého redukci kapitoly 04). Pro CS-čtoucí publikum je Podolský §6.2–§6.4 přípravně-souřadnicovým čtením; abstrakce na úroveň Lieových grup zde je modernější přeformulování spojující elementární souřadnicový aparát s Euler-Poincarého redukcí, již aplikuje kapitola 04.

---

**Oddíl 1.8 — Otevřené otázky v této kapitole** (zástupné, k vyřešení v následujících autorských průchodech):

- **OQ-1.1.** Učinit globální akci $SE(3)$ na $Q$ z §1.3.1 plně explicitní, včetně vzorce pro akci na orientacích $q_A, q_B$ parametrizovaných kvaterniony.
- **OQ-1.2.** Podat formální definici $TQ$ jako hladké variety, včetně map v okolí libovolného bodu. (Nyní popsáno neformálně.)
- **OQ-1.3.** Poskytnout odvození $\dot{\mathbf{R}}_q = \hat{\boldsymbol{\omega}}\, \mathbf{R}_q$ z prvních principů, včetně souvislosti mezi levo-invariantními a pravo-invariantními vektorovými poli na $SO(3)$. Odkaz: Marsden & Ratiu §9.

---

**Následující kapitola:** [02 — Kinetická energie na $SE(3) \times SE(3)$ v intrinzické formě](02-kineticka-energie.md) *(připravováno)*

**Konec 01-konfiguracni-varieta.md**
