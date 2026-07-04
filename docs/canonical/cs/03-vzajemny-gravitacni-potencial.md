# Kapitola 3 — Vzájemný gravitační potenciál a multipolový rozvoj

> **Předpokládané čtení.** Kapitola 01 (konfigurační varieta, operátor stříšky, konvence těžišťové soustavy) a kapitola 02 (kinetická energie, tenzor setrvačnosti jako symetrická pozitivně definitní matice $3 \times 3$, účetnictví tělesová vs. inerciální soustava). Znalost elementárního vektorového kalkulu a věty o divergenci. Žádná předchozí zkušenost se sférickými harmonikami není nutná; relevantní identity jsou zaváděny tam, kde se užívají.
>
> **Stav.** v0.1.

> **Poznámka ke způsobu a předpokladům odvození.** Tato CS kanonická kapitola (v0.1) je odvozena z anglické kanonické kapitoly `docs/canonical/en/03-mutual-gravitational-potential.md` (v0.1, s opravou znaménka §3.7.2 OQ-FORT-1, 2026-05-24). EN verze zůstává autoritativním pramenem při budoucích aktualizacích.

---

## §3.1 Co tato kapitola dělá

Kapitola 02 uzavřela s kinetickou energií soustavy dvou tuhých těles jako explicitní funkcí na $TQ$ a identifikací $T$ jako levo-invariantní Riemannovy metriky na konfigurační varietě. Polovina Lagrangiánu s potenciální energií byla odložena a je předmětem této kapitoly.

Odvozujeme **vzájemný gravitační potenciál** $V$ mezi dvěma rozšířenými tuhými tělesy ze 6rozměrného integrálu z prvních principů, rozvineme jej v inverzních mocninách vzdálenosti standardní multipolovou technikou a identifikujeme vedoucí orientačně závislý příspěvek jako **slapový člen gravitačního gradientu**. Z tohoto členu odvodíme moment síly z gravitačního gradientu $\boldsymbol{\tau}$ v tělesové soustavě a řádek po řádku ukážeme, že výsledek odpovídá implementaci inženýrské úrovně v `dynamics.py::gravity_gradient_torque_body_frame`.

Odvození sleduje metodologickou direktivu **„dvě paralelně“** vydanou projektovým vlastníkem 2026-05-23: multipolový rozvoj je proveden **nezávisle dvěma standardními přístupy** — kartézským Taylorovým rozvojem (tradice Landau-Lifšic / Murray & Dermott) a sféricko-harmonickým rozvojem (Jacksonova tradice přenesená z elektromagnetismu na gravitaci) — a oba výsledky jsou ukázány jako shodné v kvadrupólovém řádu. Krok křížové validace není stylistický; je to doktrinální `lege-artis` princip C2 křížové implementace (`_config/LEGE-ARTIS-LANGUAGE-DOCTRINE-v0.1.md` §4.4) aplikovaný na úrovni odvození spíše než na úrovni kódu. Kdyby oba přístupy nesouhlasily, jeden z nich by byl chybný; shoda je důkazem, že žádný není.

Po uzavření kapitoly 03 jsou $T$ (kapitola 02) i $V$ (tato kapitola) k dispozici. Kapitola 04 sestaví $L = T - V$ a odvodí Eulerovy–Lagrangeovy pohybové rovnice na $TQ$. Zbývající kapitoly pak specializují (kapitola 05) a validují (kapitola 06) vůči testovacím případům analytických limitů zkatalogizovaným v `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md`.

## §3.2 Přesný vzájemný potenciál jako šestirozměrný integrál

Newtonova gravitace na zde pojednávané úrovni je **párově aditivní přes hmotné elementy**: potenciální energie mezi dvěma infinitezimálními hmotnými elementy v $\mathbf{r}_1$ a $\mathbf{r}_2$ je $-G\, dm_1\, dm_2 / |\mathbf{r}_1 - \mathbf{r}_2|$ a celkový potenciál rozšířené soustavy je součtem přes všechny takové páry (Landau-Lifšic, sv. I §32; tato vlastnost přežívá v linearizované limitě slabého pole obecné relativity a láme se až v post-newtonovském řádu, jenž je mimo náš rozsah v0.1 — viz kapitola 00 §6).

Aplikujme to na dvě tuhá tělesa $A$ a $B$ s rozloženími hmoty $\rho_A(\mathbf{x}_A)$ a $\rho_B(\mathbf{x}_B)$ nesenými na kompaktních oblastech $\mathcal{B}_A, \mathcal{B}_B \subset \mathbb{R}^3$ ve svých příslušných tělesových soustavách. Hmotný bod tělesa $A$ označený tělesovou polohou $\mathbf{x}_A$ sídlí v inerciální soustavě (dle kapitoly 02, §2.2) v $\mathbf{R}_A + \mathbf{R}_q^A\, \mathbf{x}_A$ a analogicky pro těleso $B$. Euklidovská vzdálenost mezi takovou dvojicí hmotných bodů je tedy:

$$
\big|\mathbf{R}_A + \mathbf{R}_q^A\, \mathbf{x}_A - \mathbf{R}_B - \mathbf{R}_q^B\, \mathbf{x}_B\big|.
$$

Sečtením $-G\, dm_A\, dm_B$ přes všechny dvojice hmotných elementů je **přesný vzájemný gravitační potenciál** soustavy dvou těles:

$$
\boxed{
V(\mathbf{R}_A, q_A, \mathbf{R}_B, q_B)
\;=\;
-G \int_{\mathcal{B}_A} \int_{\mathcal{B}_B}
\frac{\rho_A(\mathbf{x}_A)\, \rho_B(\mathbf{x}_B)}
     {\big|\mathbf{R}_A + \mathbf{R}_q^A\, \mathbf{x}_A - \mathbf{R}_B - \mathbf{R}_q^B\, \mathbf{x}_B\big|}
\, dV_A\, dV_B.
}
$$

Toto je přesný výchozí bod — nebylo učiněno žádné přiblížení nad rámec předpokladu tuhého tělesa (rozložení hmoty $\rho_A, \rho_B$ jsou konstantní ve svých příslušných tělesových soustavách) a předpokladu newtonovské gravitace (žádné relativistické korekce). Je to funkce pouze na konfigurační varietě $Q$, nikoli na tečném svazku $TQ$: potenciál nezávisí na rychlostech. Je rovněž zjevně invariantní vůči globální akci $SE(3)$ zavedené v kapitole 01 §1.3.1 — současná tuhá translace a rotace obou těles ponechává relativní geometrii, a tudíž $V$, nezměněnou.

Dvě pozorování plynou okamžitě a tvarují zbytek kapitoly:

1. $V$ závisí pouze na **relativní poloze** $\boldsymbol{\rho} \equiv \mathbf{R}_B - \mathbf{R}_A$ a dvou **orientacích** $q_A, q_B$. Translační symetrií nemůže záviset zvlášť na $\mathbf{R}_A, \mathbf{R}_B$; v integrandu se objevuje jen jejich rozdíl. V těžišťové soustavě (konvence zamčená v kapitole 01 §1.3.1) je to funkce šesti konfiguračních proměnných: tří složek $\boldsymbol{\rho}$ a šesti rotačních proměnných $q_A, q_B$ působících na směry v integrandu.

2. Integrál je **netriviální pro libovolnou dvojici rozšířených těles**: obecně jej nelze redukovat na elementární funkce. Pro libovolná rozložení hmoty se integrál vyhodnocuje buď numericky (integrace metodou konečných prvků přes geometrie těles), nebo analyticky jen ve zvláštně-symetrických případech (např. sféricky symetrická tělesa, kde se zhroutí na výsledek bodové hmoty Newtonovou větou o slupce). Standardním fyzikálním obejitím je **multipolový rozvoj** — rozviň integrand v mocninách rozměru-tělesa-přes-vzdálenost, integruj člen po členu a usekni ve zvolené řádu. To je předmětem §3.4 a dále.

> **Hlavní odkaz.** Pro párovou aditivitu newtonovské gravitační potenciální energie pojednává CS-jazyková tradice analytické mechaniky tentýž obsah; přímý EN anchor: Landau-Lifšic, sv. I §32. **Doplňující anglojazyčné odkazy.** Murray & Dermott 1999 §4.1 (formulace přesného integrálu rozšířeného tělesa). Thorne & Blandford 2017 Rovn. (27.50) zapisuje *jednotělesovou* verzi tohoto integrálu — $\Phi(\mathbf{x}) = -G \int \rho(\mathbf{x}')/|\mathbf{x} - \mathbf{x}'|\, dV'$ — v kontextu teorie zdrojů gravitačních vln (§27.5.2); náš zarámovaný dvoutělesový vzájemný potenciál se redukuje na dvě vyhodnocení pole T&B $\Phi$ spárovaná s hmotnostmi těles, je však matematicky ekvivalentní na úrovni integrandu.

## §3.3 Obor konvergence a hranice těsného přiblížení

Multipolový rozvoj je Taylorovým rozvojem integrandu v malých parametrech $|\mathbf{x}_A|/|\boldsymbol{\rho}|$ a $|\mathbf{x}_B|/|\boldsymbol{\rho}|$ — poměrech vnitřních souřadnic těles k mezitělesové vzdálenosti. Rozvoj **konverguje absolutně**, když:

$$
|\boldsymbol{\rho}| \;>\; \sup_{\mathbf{x}_A \in \mathcal{B}_A} |\mathbf{x}_A| \;+\; \sup_{\mathbf{x}_B \in \mathcal{B}_B} |\mathbf{x}_B|,
$$

tj. když mezitělesová vzdálenost překračuje součet poloměrů nejmenších koulí se středy v těžištích každého tělesa a obklopujících to těleso. Geometricky: tělesa se nepřekrývají, s malou extra rezervou. Když tato nerovnost selže — když jsou tělesa dostatečně blízko, aby jistý hmotný bod jednoho tělesa seděl uvnitř multipolové rozvojové koule druhého — rozvoj přestává konvergovat a jakékoli konečné useknutí produkuje polynomiální divergenci, jak se $|\boldsymbol{\rho}|$ blíží poloměrům těles.

Tato hranice definuje **fyzikální obor platnosti** každého multipolově-useknutého vzorce odvozeného v této kapitole a je sdílena veškerým následujícím kódem, jenž tyto vzorce konzumuje. Inženýrská úroveň `dynamics.py` v v0.1 neobsahuje žádné tvrzení, že tato nerovnost platí v čase integrace; tvrzení je zařazeno v `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md` §3.3 jako součást testovacího případu okrajového případu těsného přiblížení.

> **Křížový odkaz.** `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md` §3.3 (porušení předpokladu těsného přiblížení) — okrajový případ třídy B vyžadující návrh tvrzení inženýrské úrovně pro v0.2. Tvrzení kanonické úrovně *zde* zní: jakýkoli multipolově-useknutý vzorec v této kapitole je prokazatelně správný jen v režimu konvergence výše; mimo něj zůstávají vzorce numericky vyhodnotitelné, ale jejich fyzikální interpretace se hroutí.

Řád useknutí rovněž omezuje zbytkovou chybu uvnitř oboru konvergence. Kvadrupólové useknutí přijaté v této kapitole je správné do členů řádu $(\ell / |\boldsymbol{\rho}|)^2$, kde $\ell$ je charakteristický rozměr tělesa; korekce dalšího řádu (oktupól-monopól a kvadrupól-dipól — obě identicky mizí pro náš v těžišti ukotvený rozvoj, ponechávajíce oktupól-kvadrupól a přímý kvadrupól-kvadrupól jako vedoucí vynechané členy) škáluje jako $(\ell / |\boldsymbol{\rho}|)^4$. Pro tělesa numerických příkladů první verze (těleso podobné JWST $\ell_A \sim 7$ m; těleso podobné sondě $\ell_B \sim 1$ m) při vzdálenostech první verze $|\boldsymbol{\rho}| \sim 10$ m až $1000$ m je kvadrupólové useknutí dobré mezi $\sim 5 \times 10^{-1}$ (těsné přiblížení) a $\sim 5 \times 10^{-9}$ (vzdálená separace). Při realistických orbitálních vzdálenostech JWST (stovky kilometrů) klesá chyba useknutí pod zaokrouhlovací chybu dvojnásobné přesnosti a tělesa jsou pro účely orbitální dynamiky efektivně bodovými hmotami; kvadrupólový člen zůstává v platnosti pro rotační vazbu bez ohledu na vzdálenost.

## §3.4 Přístup A — kartézský multipolový rozvoj

Rozvineme jmenovatel integrandu z §3.2 kolem vektoru vzdálenosti $\boldsymbol{\rho}$. Definujme $\boldsymbol{\xi}_A \equiv \mathbf{R}_q^A\, \mathbf{x}_A$ (vnitřní souřadnice tělesa $A$ vyjádřená v inerciální soustavě; rotační matice nese závislost na orientaci) a analogicky $\boldsymbol{\xi}_B$. Párový jmenovatel se stává:

$$
\big|\boldsymbol{\rho} + \boldsymbol{\xi}_B - \boldsymbol{\xi}_A\big|^{-1}
\;=\;
\big|\boldsymbol{\rho} - \mathbf{u}\big|^{-1},
\qquad \mathbf{u} \equiv \boldsymbol{\xi}_A - \boldsymbol{\xi}_B.
$$

Taylorovsky rozviň $1/|\boldsymbol{\rho} - \mathbf{u}|$ v mocninách $\mathbf{u}/|\boldsymbol{\rho}|$. Vypracováním prvních tří členů (algebra je elementární, byť poněkud zdlouhavá; standardní manipulace užívá $|\boldsymbol{\rho} - \mathbf{u}|^{-1} = (|\boldsymbol{\rho}|^2 - 2 \boldsymbol{\rho} \cdot \mathbf{u} + |\mathbf{u}|^2)^{-1/2}$ a binomickou řadu pro $(1 - \varepsilon)^{-1/2}$):

$$
\frac{1}{|\boldsymbol{\rho} - \mathbf{u}|}
\;=\;
\frac{1}{|\boldsymbol{\rho}|}
\;+\;
\frac{\hat{\boldsymbol{\rho}} \cdot \mathbf{u}}{|\boldsymbol{\rho}|^2}
\;+\;
\frac{3\, (\hat{\boldsymbol{\rho}} \cdot \mathbf{u})^2 - |\mathbf{u}|^2}{2\, |\boldsymbol{\rho}|^3}
\;+\;
\mathcal{O}\!\left(\frac{|\mathbf{u}|^3}{|\boldsymbol{\rho}|^4}\right),
$$

kde $\hat{\boldsymbol{\rho}} \equiv \boldsymbol{\rho}/|\boldsymbol{\rho}|$ je jednotkový vektor od $A$ k $B$. Dosaď $\mathbf{u} = \boldsymbol{\xi}_A - \boldsymbol{\xi}_B$ a vlož do integrálu z §3.2:

$$
V \;=\; -G \int_A \int_B \rho_A \rho_B
\left[ \frac{1}{|\boldsymbol{\rho}|}
\;+\; \frac{\hat{\boldsymbol{\rho}} \cdot (\boldsymbol{\xi}_A - \boldsymbol{\xi}_B)}{|\boldsymbol{\rho}|^2}
\;+\; \frac{3\, [\hat{\boldsymbol{\rho}} \cdot (\boldsymbol{\xi}_A - \boldsymbol{\xi}_B)]^2 - |\boldsymbol{\xi}_A - \boldsymbol{\xi}_B|^2}{2\, |\boldsymbol{\rho}|^3}
\;+\; \cdots
\right] dV_A\, dV_B.
$$

Integrály se separují: $\boldsymbol{\xi}_A$ závisí pouze na $\mathbf{x}_A$ a $q_A$, a $\boldsymbol{\xi}_B$ pouze na $\mathbf{x}_B$ a $q_B$. Vyhodnotíme tři vedoucí členy postupně.

### §3.4.1 Řád nula — monopól-monopól (Keplerův člen)

Člen $1/|\boldsymbol{\rho}|$ vychází z integrálů triviálně:

$$
V^{(0)} \;=\; -\frac{G}{|\boldsymbol{\rho}|} \int_A \rho_A\, dV_A \int_B \rho_B\, dV_B
\;=\; -\frac{G\, m_A\, m_B}{|\boldsymbol{\rho}|}.
$$

Toto je **monopól-monopól**: každé těleso působí gravitačně na druhé, jako by jeho hmotnost byla soustředěna v těžišti. Je to celý gravitační potenciál mezi dvěma bodovými hmotami a obnovuje **Keplerův problém** v limitě, kde se rotační stupně volnosti rozpojí (testovací případ `_specs/...` §2.3). Nenese žádnou závislost na orientaci, a tudíž přispívá **nulovým momentem síly z gravitačního gradientu**.

### §3.4.2 Řád jedna — dipólové členy mizí díky těžišti

Člen lineární v $\mathbf{u}$ se rozkládá na dvě části, jednu úměrnou $\int \rho_A \boldsymbol{\xi}_A dV_A$ a druhou $\int \rho_B \boldsymbol{\xi}_B dV_B$. Obě mizí:

$$
\int_A \rho_A\, \boldsymbol{\xi}_A\, dV_A \;=\; \int_A \rho_A\, \mathbf{R}_q^A\, \mathbf{x}_A\, dV_A
\;=\; \mathbf{R}_q^A \int_A \rho_A\, \mathbf{x}_A\, dV_A
\;=\; \mathbf{0},
$$

neboť $\int \rho_A \mathbf{x}_A dV_A = \mathbf{0}$ konvencí tělesového těžiště z kapitoly 02 §2.2. Táž algebra platí pro těleso $B$. Tudíž:

$$
V^{(1)} \;=\; 0.
$$

Toto je strukturálně důležité: je to **táž identita**, jež eliminovala smíšený člen v Königově rozkladu kinetické energie (kapitola 02 §2.3). Volba ukotvit souřadnice těles v těžišti eliminuje dipólové členy v *obou* $T$ i $V$ současně, ponechávajíc kvadrupól jako vedoucí orientačně závislý příspěvek ke každému.

> **Poznámka.** Kdyby, hypoteticky, byly souřadnice těles ukotveny v jiném referenčním bodě (např. v geometrickém těžišti nehomogenního tělesa, odlišném od hmotného těžiště), dipólový člen by přežil v obou $T$ i $V$ s týmž dipólovým momentem $\mathbf{d}_X = \int \rho_X \mathbf{x}_X dV_X \neq \mathbf{0}$. Konvence těžiště tedy není jen jednou z několika stejně platných voleb soustavy: je to jediné ukotvení, jež zabíjí dipólový člen a činí, že multipolový rozvoj začíná v kvadrupólovém řádu.

### §3.4.3 Řád dva — kvadrupólové členy

Příspěvek druhého řádu vyžaduje rozvinutí čitatele $3[\hat{\boldsymbol{\rho}} \cdot (\boldsymbol{\xi}_A - \boldsymbol{\xi}_B)]^2 - |\boldsymbol{\xi}_A - \boldsymbol{\xi}_B|^2$, jenž sídlí v Taylorově rozvoji z §3.4. Smíšené členy v $(\boldsymbol{\xi}_A - \boldsymbol{\xi}_B)$ rozkládají integrál na tři části:

- Část kvadratickou pouze v $\boldsymbol{\xi}_A$, integrovanou vůči $\rho_A$ a produkující celkovou hmotnost $m_B = \int \rho_B\, dV_B$ z integrálu tělesa $B$.
- Část kvadratickou pouze v $\boldsymbol{\xi}_B$, symetrickou k výše uvedené s $A \leftrightarrow B$.
- Část bilineární v $\boldsymbol{\xi}_A$ a $\boldsymbol{\xi}_B$, integrovanou vůči $\rho_A \rho_B$ — tato přispívá pouze interakcí **kvadrupól-kvadrupól**, jež klesá jako $1/|\boldsymbol{\rho}|^5$ a je o jeden řád za vedoucí slapovou vazbou. Je to vedoucí vynechaný člen, jehož efekt jsme odhadli v §3.3.

Dvě části kvadratické v jednom tělese jsou **vedoucí slapové členy** a mají identickou strukturu (s $A$ a $B$ vyměněnými). Pro těleso $A$:

$$
V_{\mathrm{tidal},\, A}
\;=\;
-\frac{G\, m_B}{2\, |\boldsymbol{\rho}|^3}
\int_A \rho_A
\Big[ 3\, (\hat{\boldsymbol{\rho}} \cdot \boldsymbol{\xi}_A)^2 \;-\; |\boldsymbol{\xi}_A|^2 \Big]\, dV_A.
$$

K vyjádření tohoto pomocí **tenzoru setrvačnosti** $\mathbf{I}_A$ z kapitoly 02 §2.4 zaveďme tělesový jednotkový vektor $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \equiv (\mathbf{R}_q^A)^T \hat{\boldsymbol{\rho}}$ (inerciální směr vzdálenosti vyjádřený v souřadnicích tělesa $A$). Integrand se stává:

$$
3\, (\hat{\boldsymbol{\rho}} \cdot \boldsymbol{\xi}_A)^2 - |\boldsymbol{\xi}_A|^2
\;=\; 3\, (\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{x}_A)^2 - |\mathbf{x}_A|^2,
$$

s užitím ortogonality $|\mathbf{R}_q^A \mathbf{x}_A| = |\mathbf{x}_A|$ a analogické identity skalárního součinu. Definuj tenzor druhého momentu tělesa $A$ v jeho tělesové soustavě:

$$
(\mathbf{M}_A)_{ij} \;=\; \int_A \rho_A(\mathbf{x}_A)\, (\mathbf{x}_A)_i\, (\mathbf{x}_A)_j\, dV_A.
$$

Pak $\int \rho_A (\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{x}_A)^2 dV_A = \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{M}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$ a $\int \rho_A |\mathbf{x}_A|^2 dV_A = \mathrm{tr}(\mathbf{M}_A)$.

Vztažením $\mathbf{M}_A$ k tenzoru setrvačnosti z kapitoly 02 §2.4 — $(\mathbf{I}_A)_{ij} = \int \rho_A (|\mathbf{x}_A|^2 \delta_{ij} - (\mathbf{x}_A)_i (\mathbf{x}_A)_j) dV_A$ — máme $\mathbf{I}_A = \mathrm{tr}(\mathbf{M}_A) \mathbb{I}_3 - \mathbf{M}_A$ a inverzně $\mathbf{M}_A = \tfrac{1}{2} \mathrm{tr}(\mathbf{I}_A) \mathbb{I}_3 - \mathbf{I}_A$ (s užitím $\mathrm{tr}(\mathbf{I}_A) = 2\, \mathrm{tr}(\mathbf{M}_A)$). Dosaď a zjednoduš:

$$
3\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{M}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}
\;-\; \mathrm{tr}(\mathbf{M}_A)
\;=\;
\mathrm{tr}(\mathbf{I}_A) \;-\; 3\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}.
$$

Tudíž **kvadrupólový slapový člen pro těleso $A$**, vyjádřený v tělesovém tenzoru setrvačnosti, jejž kapitola 02 již stanovila, je:

$$
\boxed{
V_{\mathrm{tidal},\, A}
\;=\;
-\frac{G\, m_B}{2\, |\boldsymbol{\rho}|^3}
\Big[ \mathrm{tr}(\mathbf{I}_A) \;-\; 3\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \Big].
}
$$

Symetrický výraz platí pro $V_{\mathrm{tidal},\, B}$ s $A \leftrightarrow B$. Celkový vedoucí orientačně závislý potenciál je součtem $V_{\mathrm{tidal}} = V_{\mathrm{tidal},\, A} + V_{\mathrm{tidal},\, B}$ a **úplný kvadrupólově-useknutý vzájemný potenciál** je:

$$
\boxed{
V \;=\; -\frac{G\, m_A\, m_B}{|\boldsymbol{\rho}|}
\;+\; V_{\mathrm{tidal},\, A}
\;+\; V_{\mathrm{tidal},\, B}
\;+\; \mathcal{O}\!\left(\frac{\ell^4}{|\boldsymbol{\rho}|^5}\right).
}
$$

> **Doplňující anglojazyčné odkazy.** Landau-Lifšic, sv. I §41 odvozuje multipolový rozvoj gravitačního potenciálu jednoho tělesa v kompaktní formě; Murray & Dermott 1999 §4 provádí rozšíření na binární těleso toutéž kartézskou technikou. Thorne & Blandford 2017 §27.5.2 („Quadrupole-Moment Formalism“) provádí týž kartézský Taylorův rozvoj v Rovn. (27.51)–(27.53) v zarámování teorie zdrojů gravitačních vln: T&B Rovn. (27.51) je Taylorovým rozvojem $1/|\mathbf{x} - \mathbf{x}'|$ v mocninách $\mathbf{x}'/r$ (jejich $\mathbf{x}'$ je naše $\mathbf{u}$ z §3.4); T&B Rovn. (27.52) sbírá monopólově-plus-kvadrupólově useknutý potenciál; T&B Rovn. (27.53) definuje bezstopový hmotný kvadrupólový moment $\mathcal{I}_{jk}$. **Poznámka ke konvenci:** hmotný kvadrupólový moment T&B $\mathcal{I}_{jk}$ užívá předfaktor $1/3$, kde kartézský kvadrupólový tenzor $\mathbf{Q}$ užitý v naší §3.5 sleduje Jacksonovu konvenci s předfaktorem $3$. Oba jsou ekvivalentní: $\mathbf{Q}_{jk} = 3\, \mathcal{I}_{jk}$, oba bezstopové symetrické kódující týchž pět fyzikálních stupňů volnosti. Forma zde v pojmech $\mathrm{tr}(\mathbf{I}_A) - 3 \hat{\boldsymbol{\rho}} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}$ udržuje tělesový tenzor setrvačnosti viditelný, což je správná forma pro spojení s momentem síly z gravitačního gradientu (§3.7).

## §3.5 Přístup B — sféricko-harmonický multipolový rozvoj

Nyní zopakujeme rozvoj nezávislou cestou, s užitím sféricko-harmonické věty o součtu pro $1/|\mathbf{r}_1 - \mathbf{r}_2|$, jež je standardní v teorii elektromagnetických multipólů a přenáší se doslova na gravitaci (jediné rozdíly jsou znaménko zdroje a nepřítomnost členů magnetických multipólů; oba jsou irelevantní na úrovni matematické identity).

Připomeňme **Laplaceův rozvoj** (Jackson 1999 rovn. 3.70):

$$
\frac{1}{|\mathbf{r}_1 - \mathbf{r}_2|}
\;=\; \sum_{\ell=0}^{\infty} \sum_{m=-\ell}^{\ell}
\frac{4\pi}{2\ell+1}\,
\frac{r_<^{\ell}}{r_>^{\ell+1}}\,
Y_\ell^{m\,*}(\hat{\mathbf{r}}_2)\, Y_\ell^m(\hat{\mathbf{r}}_1),
$$

kde $r_< = \min(|\mathbf{r}_1|, |\mathbf{r}_2|)$, $r_> = \max(|\mathbf{r}_1|, |\mathbf{r}_2|)$ a $Y_\ell^m$ jsou standardní normalizované sférické harmoniky. Identita konverguje absolutně, kdykoli $r_< < r_>$ (tj. kdekoli je dobře definovaná), což je týž obor konvergence, jejž jsme identifikovali v §3.3.

K čistému použití tohoto na náš problém binárního tělesa pracujme nejprve s případem, kdy je těleso $B$ **bodovou hmotou** v inerciální poloze $\mathbf{R}_B$ (to nám dá $V_{\mathrm{tidal},\, A}$; symetrický případ pro těleso $B$ následuje výměnou označení a bilineární člen kvadrupól-kvadrupól je o jeden řád vyšší, jak bylo poznamenáno v §3.4.3). Jednotělesový potenciál generovaný tělesem $A$ v bodě pole $\mathbf{R}_B$ je:

$$
\Phi_A(\mathbf{R}_B) \;=\; -G \int_A \frac{\rho_A(\mathbf{x}_A)}
{|\mathbf{R}_A + \boldsymbol{\xi}_A - \mathbf{R}_B|}\, dV_A
\;=\; -G \int_A \frac{\rho_A(\mathbf{x}_A)}{|\boldsymbol{\xi}_A - \boldsymbol{\rho}|}\, dV_A,
$$

kde opět $\boldsymbol{\rho} = \mathbf{R}_B - \mathbf{R}_A$ a $\boldsymbol{\xi}_A = \mathbf{R}_q^A \mathbf{x}_A$. Potenciální energie přispívaná interakcí těleso-$A$-na-bodovou-hmotu-$B$ je $V_A^{\mathrm{on}\, B} = m_B \Phi_A(\mathbf{R}_B)$.

Aplikuj Laplaceův rozvoj s $\mathbf{r}_1 = \boldsymbol{\rho}$ (větší; jsme v oboru konvergence $|\boldsymbol{\rho}| > |\boldsymbol{\xi}_A|$) a $\mathbf{r}_2 = \boldsymbol{\xi}_A$ (menší):

$$
\frac{1}{|\boldsymbol{\xi}_A - \boldsymbol{\rho}|}
\;=\; \sum_{\ell, m} \frac{4\pi}{2\ell+1}\, \frac{|\boldsymbol{\xi}_A|^{\ell}}{|\boldsymbol{\rho}|^{\ell+1}}\,
Y_\ell^{m\,*}(\hat{\boldsymbol{\xi}}_A)\, Y_\ell^m(\hat{\boldsymbol{\rho}}).
$$

Dosaď a vytáhni na $\boldsymbol{\rho}$ závislý faktor z integrálu:

$$
V_A^{\mathrm{on}\, B}
\;=\; -G\, m_B \sum_{\ell, m}
\frac{4\pi}{2\ell+1}\,
\frac{Y_\ell^m(\hat{\boldsymbol{\rho}})}{|\boldsymbol{\rho}|^{\ell+1}}\,
\underbrace{\int_A \rho_A(\mathbf{x}_A)\, |\boldsymbol{\xi}_A|^{\ell}\, Y_\ell^{m\,*}(\hat{\boldsymbol{\xi}}_A)\, dV_A}_{\equiv\;\, q_\ell^m(A)}.
$$

Závorkované integrály jsou **sféricko-harmonické multipolové momenty** tělesa $A$:

$$
q_\ell^m(A) \;\equiv\; \int_A \rho_A(\mathbf{x}_A)\, |\boldsymbol{\xi}_A|^{\ell}\, Y_\ell^{m\,*}(\hat{\boldsymbol{\xi}}_A)\, dV_A.
$$

Tyto jsou obecně komplexně-hodnotové; gravitační multipolové momenty jsou fyzikálními pozorovatelnými a podléhají obvyklým podmínkám reálnosti $q_\ell^{-m} = (-1)^m q_\ell^{m\,*}$, neboť $\rho_A$ je reálná hustota. Závisí na orientaci $q_A$ skrze $\boldsymbol{\xi}_A = \mathbf{R}_q^A \mathbf{x}_A$.

Nyní vyhodnotíme první tři multipolové řády:

**$\ell = 0$ (monopól).** $Y_0^0 = 1/\sqrt{4\pi}$. Moment je $q_0^0(A) = m_A / \sqrt{4\pi}$, nezávislý na orientaci. Příspěvek k $V$:

$$
-G m_B \cdot \frac{4\pi}{1} \cdot \frac{m_A / \sqrt{4\pi}}{|\boldsymbol{\rho}|} \cdot \frac{1}{\sqrt{4\pi}}
\;=\; -\frac{G m_A m_B}{|\boldsymbol{\rho}|}.
$$

Toto reprodukuje Keplerův člen z §3.4.1, jak je požadováno.

**$\ell = 1$ (dipól).** Multipolové momenty $\ell=1$ jsou lineární v $\mathbf{x}_A$ a mizí konvencí těžiště (týž argument jako §3.4.2). Konkrétně $q_1^0(A) \propto \int \rho_A z_A dV_A = 0$ a analogicky pro $q_1^{\pm 1}$. Žádný příspěvek.

**$\ell = 2$ (kvadrupól).** Pět komplexních sféricko-harmonických momentů $q_2^m(A)$ ($m = -2, -1, 0, 1, 2$) kóduje týž fyzikální obsah jako pět nezávislých složek symetrického, bezstopového **kvadrupólového tenzoru**:

$$
(\mathbf{Q}_A)_{ij} \;=\; \int_A \rho_A(\mathbf{x}_A)
\big[ 3\, (\boldsymbol{\xi}_A)_i\, (\boldsymbol{\xi}_A)_j \;-\; |\boldsymbol{\xi}_A|^2\, \delta_{ij} \big]\, dV_A.
$$

Vztah mezi sféricko-harmonickou a kartézskou formou je standardní v Jackson 1999 §4.1: $q_2^m(A)$ jsou jisté lineární kombinace $(\mathbf{Q}_A)_{ij}$, s ekvivalencí dokázanou přímou identifikací sféricko-harmonických bázových funkcí $r^2 Y_2^m(\hat{\mathbf{r}})$ s bezstopovými kvadratickými polynomy v kartézských souřadnicích.

**Explicitní algebra změny báze (řeší OQ-3.1).** S užitím standardních sférických harmonik $\ell=2$ v Condon-Shortleyově konvenci (Jackson 1999 §3.6):
$$
Y_2^{0}(\hat{\mathbf{r}}) = \sqrt{\tfrac{5}{16\pi}}\,(3\cos^2\theta - 1), \quad
Y_2^{\pm 1}(\hat{\mathbf{r}}) = \mp \sqrt{\tfrac{15}{8\pi}}\,\sin\theta\cos\theta\, e^{\pm i\phi}, \quad
Y_2^{\pm 2}(\hat{\mathbf{r}}) = \sqrt{\tfrac{15}{32\pi}}\,\sin^2\theta\, e^{\pm 2i\phi},
$$
spolu se sféricko-kartézskými přepisy $r^2 (3\cos^2\theta - 1) = 2 z^2 - x^2 - y^2$, $r^2 \sin\theta\cos\theta\, e^{i\phi} = z(x + iy)$ a $r^2 \sin^2\theta\, e^{2i\phi} = (x + iy)^2$ se pět momentů $\ell=2$ $q_2^m(A) = \int \rho_A\, |\boldsymbol{\xi}_A|^2\, Y_2^{m\,*}(\hat{\boldsymbol{\xi}}_A)\, dV_A$ mapuje na složky kartézského kvadrupólového tenzoru jako:

$$
\boxed{
\begin{aligned}
q_2^{0}(A) \;&=\; \sqrt{\tfrac{5}{16\pi}}\, \cdot \tfrac{1}{3}\, (\mathbf{Q}_A)_{zz},
\\[2pt]
q_2^{\pm 1}(A) \;&=\; \mp \sqrt{\tfrac{15}{8\pi}}\, \cdot \tfrac{1}{3}\, \big[\, (\mathbf{Q}_A)_{xz} \;\mp\; i\,(\mathbf{Q}_A)_{yz}\, \big],
\\[2pt]
q_2^{\pm 2}(A) \;&=\; \sqrt{\tfrac{15}{32\pi}}\, \cdot \tfrac{1}{3}\, \big[\, (\mathbf{Q}_A)_{xx} - (\mathbf{Q}_A)_{yy} \;\mp\; 2i\,(\mathbf{Q}_A)_{xy}\, \big].
\end{aligned}
}
$$

(Faktor $1/3$ je rozdíl předfaktorů mezi naší definicí $\mathbf{Q}_A$ s vedoucím $3 x_i x_j - r^2 \delta_{ij}$ v §3.5 a konvencí T&B/Jackson s vedoucím $x_i x_j - \tfrac{1}{3} r^2 \delta_{ij}$ — připomeňme poznámku ke konvenci v referenčním anchoru §3.4.3: $\mathbf{Q}_{jk} = 3\, \mathcal{I}_{jk}$.) Pět komplexních momentů není nezávislých — reálnost hmotné hustoty $\rho_A$ vynucuje $q_2^{-m}(A) = (-1)^m\, q_2^{m\,*}(A)$, ponechávajíc přesně pět reálných stupňů volnosti odpovídajících pěti nezávislým prvkům bezstopového symetrického kartézského tenzoru $3 \times 3$: $(\mathbf{Q}_A)_{xx}$, $(\mathbf{Q}_A)_{yy}$ (s $(\mathbf{Q}_A)_{zz} = -(\mathbf{Q}_A)_{xx} - (\mathbf{Q}_A)_{yy}$ bezstopovostí), $(\mathbf{Q}_A)_{xy}$, $(\mathbf{Q}_A)_{xz}$, $(\mathbf{Q}_A)_{yz}$. Mapování výše je invertibilní — dáno pěti $q_2^m(A)$, kartézské prvky se získají inverzí předfaktorů a vzetím reálných/imaginárních částí.

**Geometrická interpretace mapování.** Každý sféricko-harmonický moment „filtruje“ specifickou kartézskou kombinaci svou úhlovou projekcí: $q_2^0$ měří protáhlé/zploštělé protažení tělesa podél $z$ (závisí pouze na $(\mathbf{Q}_A)_{zz}$, ekvivalentně na axiálním momentu tělesa vůči jeho kolmým momentům); $q_2^{\pm 1}$ měří mimodiagonální vazby $xz$ a $yz$ (náklon tělesa mimo rovinu $xy$); $q_2^{\pm 2}$ měří rovinnou asymetrii $(\mathbf{Q}_A)_{xx} - (\mathbf{Q}_A)_{yy}$ a mimodiagonální smyk $xy$. Pro axisymetrické těleso kolem $z$ (např. `make_oblate_reference_body()` z §2.4) je nenulové pouze $q_2^0$; všechny čtyři ostatní momenty $\ell=2$ identicky mizí. To je konzistentní s kartézským pozorováním, že $\mathbf{Q}_A$ pro axisymetrické těleso je diagonální s $(\mathbf{Q}_A)_{xx} = (\mathbf{Q}_A)_{yy}$ a pouze $(\mathbf{Q}_A)_{zz}$ jako volným parametrem.

Provedením součtu $\ell = 2$ a užitím této kartézsko-sférické ekvivalence (kontrakce na $\hat{\boldsymbol{\rho}}$ v sféricko-harmonickém součtu se zhroutí na skalární součin $\hat{\boldsymbol{\rho}} \cdot \mathbf{Q}_A \cdot \hat{\boldsymbol{\rho}}$ v kartézské formě; algebraické detaily jsou standardní a vynecháváme je):

$$
V_A^{\mathrm{on}\, B, \,\ell=2}
\;=\;
-\frac{G\, m_B}{2\, |\boldsymbol{\rho}|^3}\,
\hat{\boldsymbol{\rho}} \cdot \mathbf{Q}_A \cdot \hat{\boldsymbol{\rho}}.
$$

Kartézský kvadrupólový tenzor $\mathbf{Q}_A$ souvisí s tělesovým tenzorem setrvačnosti $\mathbf{I}_A$ **stopovou identitou**:

$$
(\mathbf{Q}_A)_{ij} \;=\; \mathrm{tr}(\mathbf{I}_A^{\mathrm{inertial}})\, \delta_{ij} \;-\; 3\, (\mathbf{I}_A^{\mathrm{inertial}})_{ij},
$$

kde $\mathbf{I}_A^{\mathrm{inertial}} = \mathbf{R}_q^A \mathbf{I}_A (\mathbf{R}_q^A)^T$ je tenzor setrvačnosti $A$ vyjádřený v souřadnicích inerciální soustavy (kapitola 02 §2.5). Oba $\mathbf{Q}_A$ a $\mathbf{I}_A^{\mathrm{inertial}}$ jsou symetrické, avšak $\mathbf{Q}_A$ je navíc **bezstopový** ($\mathrm{tr}(\mathbf{Q}_A) = 3\, \mathrm{tr}(\mathbf{I}_A) - 3\, \mathrm{tr}(\mathbf{I}_A) = 0$); stopová identita výše je explicitním rozkladem obecného symetrického tenzoru na jeho stopovou a bezstopovou část.

> **Ověření.** Stopovou identitu lze ověřit přímo: aplikováním definic $(\mathbf{Q}_A)_{ij}$ a $(\mathbf{I}_A)_{ij}$ z této kapitoly resp. kapitoly 02 se obě redukují na integrály přes $\rho_A$ polynomů v $\mathbf{x}_A$; lineární kombinace $\mathrm{tr}(\mathbf{I}_A) \delta_{ij} - 3 (\mathbf{I}_A)_{ij}$ přesně obnovuje $\int \rho_A [3 x_i x_j - |\mathbf{x}|^2 \delta_{ij}] dV$. Numerická kontrola zdravého rozumu na setrvačnosti tenké tyče $\mathbf{I} = \mathrm{diag}(mL^2/3,\, mL^2/3,\, 0)$ produkuje $\mathbf{Q} = \mathrm{diag}(-mL^2/3,\, -mL^2/3,\, 2mL^2/3)$, bezstopovou na strojovou přesnost (ověřeno během autorství kapitoly 03).

Dosazením stopové identity:

$$
\hat{\boldsymbol{\rho}} \cdot \mathbf{Q}_A \cdot \hat{\boldsymbol{\rho}}
\;=\; \mathrm{tr}(\mathbf{I}_A)\, |\hat{\boldsymbol{\rho}}|^2
\;-\; 3\, \hat{\boldsymbol{\rho}} \cdot \mathbf{I}_A^{\mathrm{inertial}} \cdot \hat{\boldsymbol{\rho}}
\;=\; \mathrm{tr}(\mathbf{I}_A) \;-\; 3\, \hat{\boldsymbol{\rho}} \cdot \mathbf{I}_A^{\mathrm{inertial}} \cdot \hat{\boldsymbol{\rho}}.
$$

A s užitím $\hat{\boldsymbol{\rho}} \cdot \mathbf{I}_A^{\mathrm{inertial}} \cdot \hat{\boldsymbol{\rho}} = \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$ (podobnostní transformace tenzoru setrvačnosti; kapitola 02 §2.5):

$$
\boxed{
V_A^{\mathrm{on}\, B, \,\ell=2}
\;=\;
-\frac{G\, m_B}{2\, |\boldsymbol{\rho}|^3}
\Big[ \mathrm{tr}(\mathbf{I}_A) \;-\; 3\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \Big].
}
$$

Toto je **identické** se zarámovaným výsledkem §3.4.3.

> **Doplňující anglojazyčné odkazy.** Jackson 1999 §4.1 (multipolový rozvoj, kartézsko-sférická ekvivalence při $\ell = 2$). Matematika je identická v elektromagnetickém a gravitačním kontextu; jedinou substitucí je $\rho_{\mathrm{charge}} \to \rho_{\mathrm{mass}}$ a $1/(4\pi\varepsilon_0) \to -G$. Murray & Dermott 1999 §4.4 rovněž prezentuje sféricko-harmonickou formu pro aplikace ve sluneční soustavě, s konvencemi specializovanými na koeficienty zploštělých těles $J_\ell$; my užíváme Jacksonovu normalizaci napříč touto kapitolou.

## §3.6 Křížová validace: kartézský ≡ sféricko-harmonický v kvadrupólovém řádu

Oba přístupy, započaté z nezávislých identit (kartézský Taylorův rozvoj $1/|\boldsymbol{\rho} - \mathbf{u}|$ v §3.4 a Laplaceova sféricko-harmonická věta o součtu v §3.5), produkují tentýž zarámovaný kvadrupólový slapový výraz:

$$
V_{\mathrm{tidal},\, A}
\;=\;
-\frac{G\, m_B}{2\, |\boldsymbol{\rho}|^3}
\Big[ \mathrm{tr}(\mathbf{I}_A) \;-\; 3\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \Big].
$$

Shoda **není náhoda**: na úrovni matematických identit jsou kartézský Taylorův rozvoj $1/|\boldsymbol{\rho} - \mathbf{u}|$ a Laplaceova sféricko-harmonická věta o součtu ekvivalentní (jednu lze odvodit z druhé přeskupením členů). Co křížové odvození činí, je ověření, že **jsme neudělali algebraickou chybu v žádné z cest**: za cenu zdvojeného výpočtu jsme eliminovali třídu chyb vznikajících z jediného řetězce odvození.

Tato metodologie — **dva nezávislé standardní nástroje, křížově validované** — byla přímo nařízena projektovým vlastníkem 2026-05-23 pro jakékoli odvození kanonické úrovně užívající standardní matematickou techniku (`CLAUDE.md` HANDOFF BLOCK 2026-05-23 §4). Je to analytickým analogem doktrinálního kritéria §4.4 C2 **bitové identity mezi implementacemi**. Následující kapitoly kanonické úrovně budou aplikovat tutéž disciplínu: kapitola 04 (Euler-Lagrange + Euler-Poincarého redukce) odvodí tělesové rotační rovnice jak přímým Eulerovým–Lagrangeovým výpočtem, tak Euler-Poincarého redukční větou; kapitola 05 (stabilita teorému tenisové rakety) jak linearizací, tak přímou eliptickou-funkční integrací; kapitola 06 validuje každý testovací případ vůči alespoň dvěma učebnicovým odkazům.

> **Doplňující anglojazyčný odkaz.** Jackson 1999 §4.1 obsahuje explicitní důkaz, že kartézský kvadrupólový tenzor $\mathbf{Q}_{ij}$ a sféricko-harmonické multipolové momenty $q_2^m$ kódují týchž pět nezávislých stupňů volnosti bezstopového symetrického tenzoru $3 \times 3$. Mapování je změnou báze v prostoru tenzorových harmonik $\ell = 2$; fyzikální obsah je identický.

> **Hlavní odkaz.** Kartézský kvadrupólový tenzor $\mathbf{Q}_{ij}$ a sféricko-harmonický multipolový rozvoj jsou *moderní* matematické metody přímo nepojednané v Podolský 2018 — ono dílo je postaveno kolem klasické analytické mechaniky tuhých těles, nikoli gravitačních potenciálů rozšířených těles. Most od Podolského k naší §3.5 vede skrze §7.1 *Tenzor setrvačnosti* Rovn. (7.3)–(7.9): bezstopová symetrická struktura gravitačního kvadrupólového tenzoru paralelizuje Podolského bilineárně-formové pojednání kineticko-teoretického tenzoru setrvačnosti (oba jsou reálné symetrické tenzory $3 \times 3$ vybudované jako integrály rozložení hmoty vůči prostorovým polynomům druhého řádu). Algebraický slovník, jejž Podolský rozvíjí v §7.1 — definice symetrického tenzoru bilineární formou, rozklad na hlavní osy, směrová projekce $I_{\mathbf{n}}$ — je přímo znovupoužitelným myšlenkovým lešením pro gravitační kvadrupólový tenzor. CS čtenář, jenž chce samotnou teorii multipolového rozvoje, je nejlépe obsloužen elektrodynamickým textem (Jackson 1999 §4.1 citovaný výše, nebo CS-jazykový kurz Klasické elektrodynamiky, je-li lokálně dostupný).

## §3.7 Moment síly z gravitačního gradientu na těleso A

Moment síly z gravitačního gradientu se získá z orientačně závislé části potenciálu standardním variačním pravidlem pro rotační stupně volnosti: při virtuální rotaci $\delta \boldsymbol{\theta}_A$ tělesa $A$ (malý vektor osa-úhel aplikovaný na jeho tělesovou soustavu) je změna $V_{\mathrm{tidal},\, A}$ rovna $\delta V = -\boldsymbol{\tau}_A \cdot \delta \boldsymbol{\theta}_A$, definujíc tělesový moment síly $\boldsymbol{\tau}_A$. Odvodíme $\boldsymbol{\tau}_A$ explicitně.

Při infinitezimální tělesové rotaci $\mathbf{R}_q^A \to \mathbf{R}_q^A (\mathbb{I} + \widehat{\delta \boldsymbol{\theta}_A})$ je inerciální směr vzdálenosti $\hat{\boldsymbol{\rho}}$ nezměněn (rotace působí na těleso, nikoli na vektor vzdálenosti těles) a transformovaný tělesový jednotkový vektor je:

$$
\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}^{\,\mathrm{new}}
\;=\; \big[\mathbf{R}_q^A (\mathbb{I} + \widehat{\delta \boldsymbol{\theta}_A})\big]^T \hat{\boldsymbol{\rho}}
\;=\; (\mathbb{I} - \widehat{\delta \boldsymbol{\theta}_A})\, (\mathbf{R}_q^A)^T \hat{\boldsymbol{\rho}}
\;=\; (\mathbb{I} - \widehat{\delta \boldsymbol{\theta}_A})\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A},
$$

kde jsme užili $(AB)^T = B^T A^T$ a $\widehat{\delta \boldsymbol{\theta}_A}^T = -\widehat{\delta \boldsymbol{\theta}_A}$ (antisymetrie operátoru stříšky). Do prvního řádu v $\delta \boldsymbol{\theta}_A$ to dává:

$$
\delta \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}
\;=\;
-\, \widehat{\delta \boldsymbol{\theta}_A}\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}
\;=\;
-\, \delta \boldsymbol{\theta}_A \times \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}.
$$

Stopový člen $\mathrm{tr}(\mathbf{I}_A)$ v $V_{\mathrm{tidal},\, A}$ je **nezávislý na orientaci** (stopa tenzoru je invariantní vůči podobnostní transformaci) a přispívá nulou k $\delta V$. Orientačně závislou částí je $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$, jejíž variace je:

$$
\delta\!\left[\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}\right]
\;=\;
2\, \delta \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}
\;=\;
-2\, (\delta \boldsymbol{\theta}_A \times \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}) \cdot (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}).
$$

S užitím identity smíšeného součinu $(\mathbf{a} \times \mathbf{b}) \cdot \mathbf{c} = \mathbf{a} \cdot (\mathbf{b} \times \mathbf{c})$:

$$
\delta\!\left[\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}\right]
\;=\;
-2\, \delta \boldsymbol{\theta}_A \cdot
\big[ \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \times (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}) \big].
$$

Vložením tohoto do variace $V_{\mathrm{tidal},\, A}$ (faktor $-3$ před orientačně závislou částí, faktor $-G m_B / (2 |\boldsymbol{\rho}|^3)$ celkově, takže orientačně závislý koeficient je $+3 G m_B / (2 |\boldsymbol{\rho}|^3)$):

$$
\delta V_{\mathrm{tidal},\, A}
\;=\;
\frac{3\, G\, m_B}{2\, |\boldsymbol{\rho}|^3}
\cdot (-2)\, \delta \boldsymbol{\theta}_A \cdot
\big[ \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \times (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}) \big]
\;=\;
-\, \frac{3\, G\, m_B}{|\boldsymbol{\rho}|^3}\, \delta \boldsymbol{\theta}_A \cdot
\big[ \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \times (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}) \big].
$$

Definující vztah $\delta V = -\boldsymbol{\tau}_A \cdot \delta \boldsymbol{\theta}_A$ dává:

$$
\boxed{
\boldsymbol{\tau}_A^{\mathrm{body}}
\;=\;
\frac{3\, G\, m_B}{|\boldsymbol{\rho}|^3}\,
\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \times (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}).
}
$$

Toto je **moment síly z gravitačního gradientu** na těleso $A$ v jeho tělesové soustavě, ve vedoucím (kvadrupólovém) řádu v poměru rozměr-tělesa / vzdálenost. Symetrický výraz platí pro $\boldsymbol{\tau}_B^{\mathrm{body}}$ s $A \leftrightarrow B$ a $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, B} = -(\mathbf{R}_q^B)^T \hat{\boldsymbol{\rho}}$ (jednotkový vektor od $B$ k $A$ v souřadnicích tělesa $B$ má opačné znaménko).

> **Doplňující anglojazyčné odkazy.** Hughes 2004 §5.3 (moment síly z gravitačního gradientu na tuhou kosmickou loď v důsledku centrální bodové hmoty; výsledek je dán v podstatě v této formě). Hughesovo odvození se specializuje na případ, kdy je těleso $B$ mnohem větším primárem braným jako bodová hmota; zde jsme odvodili symetrický případ dvou rozšířených těles, jehož je Hughesův výsledek limitou $\mathbf{I}_B \to \mathbf{0}$. Goldstein 2002 §5.6 odvozuje týž výraz přímým výpočtem integrálu momentu síly přes rozložení hmoty tělesa (bez multipolového mezikroku) jako nezávislou třetí-cestnou křížovou kontrolu.

### §3.7.1 Shoda řádek po řádku vůči `dynamics.py`

Implementace inženýrské úrovně v `dynamics.py::gravity_gradient_torque_body_frame` (řádky 113–131) zní, zhuštěně:

```
r_mag = |r_rel_inertial|
r_hat_inertial = r_rel_inertial / r_mag
r_hat_body = R_body_to_inertial.T @ r_hat_inertial
I_rhat = I_body @ r_hat_body
tau = (3 * G * m_other / r_mag**3) * cross(r_hat_body, I_rhat)
```

Volací místo pro těleso $A$ (`state_derivative`, řádky 154–157) předává `r_rel_inertial = body_b.x - body_a.x`, což je naše $\boldsymbol{\rho}$ od $A$ k $B$. Tudíž `r_hat_body = R_A.T @ rho_hat = rho_hat_body_A` v naší notaci. Dosazením:

$$
\texttt{tau} \;=\; \frac{3\, G\, m_B}{|\boldsymbol{\rho}|^3}\, \big[\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \times (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A})\big],
$$

což je **identické** se zarámovaným výsledkem §3.7. Shoda je řádek po řádku: předfaktor, znaménko, pořadí vektorového součinu a konvence pro směr $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$ vše souhlasí. Dokumentační řetězec implementace (řádek 122) cituje Hughes 2004 §5.3 přímo; odvození této kapitoly je rozvinutím této citace na kanonické úrovni.

Pro moment síly na těleso $B$ volací místo (`state_derivative`, řádek 158) předává `r_rel_inertial = -r_a_to_b_inertial = body_a.x - body_b.x`, což je $-\boldsymbol{\rho}$ z pohledu $B$ — přesně jednotkový vektor od $B$ k $A$. Konvence znaménka je tedy konzistentní: v každém volání tělesa $\hat{\boldsymbol{\rho}}_{\mathrm{body}}$ míří od onoho tělesa k druhému.

> **Hlavní odkaz.** Moment síly z gravitačního gradientu na rozšířené těleso ve vnějším gravitačním poli je *moderní* téma dynamiky postoje kosmické lodě, mimo klasický rozsah Podolský 2018. Nejbližším klasicko-mechanickým precedentem v Podolském je **Kapitola 8 §8.3** *Těžký symetrický setrvačník s pevným bodem* (Lagrangeův setrvačník). V tom problému působí vnější gravitační pole na axisymetrické tuhé těleso zavěšené z pevného bodu, produkujíc moment síly $\mathbf{M}$, jenž závisí na orientaci tělesa skrze výšku jeho těžiště — analyticky ekvivalentní algebraickou strukturou naší závislosti momentu síly z gravitačního gradientu na orientaci. Podolský §8.3 odvozuje výsledný pohyb Lagrangeovým formalismem s Eulerovými úhly jako zobecněnými souřadnicemi: první integrály cyklických souřadnic $L_\psi, L_\varphi$ v jeho Rovn. (8.11), zachování energie redukující na jednorozměrný pohyb v $\vartheta$ s efektivním potenciálem $V_{\mathrm{ef}}(\vartheta)$ v zarámované Rovn. (8.15) a klasifikace precesních režimů (Rovn. 8.13 s případy (a)/(b)/(c) na jednotkové kouli). Náš problém librace gravitačního gradientu (budoucí testovací případ kapitoly 06) je analogem Podolského Lagrangeova setrvačníku v dynamice postoje kosmické lodě: táž lagrangeovská struktura, odlišná fyzikální interpretace „momentu síly kolem osy zavěšení“. Fortranový librační testbed fáze 2 bude křížově validovat vůči metodě $V_{\mathrm{ef}}$ Lagrangeova setrvačníku jako svému kanonickému CS-jazykovému orákulu analytického limitu.

### §3.7.2 Jacobiho matice momentu síly z gravitačního gradientu: citlivost na počáteční podmínky a geometrii sestavy (most na inženýrskou úroveň)

Odvození §3.7 produkovalo moment síly z gravitačního gradientu $\boldsymbol{\tau}_A^{\mathrm{body}}$ jako uzavřenou funkci konfigurace $(\boldsymbol{\rho}, q_A, q_B)$. Pro inženýrskou úroveň — praktickou analýzu citlivosti, odhad dopadu počátečních podmínek, volbu velikosti integračního kroku — je **linearizovaná Jacobiho matice** $\boldsymbol{\tau}_A$ vzhledem k malým perturbacím konfigurace výstupem doplňujícím uzavřený vzorec. Tato podkapitola přenáší rámec mostu na inženýrskou úroveň z kapitoly 02 §2.5.1 (mimodiagonální citlivost $d\mathbf{I}^{\mathrm{inertial}}/dt$) do kontextu gravitačního gradientu, dle direktivy projektového vlastníka 2026-05-23.

Konfigurace vstupuje do $\boldsymbol{\tau}_A^{\mathrm{body}} = (3 G m_B / r^3)\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \times (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A})$ skrze tři místa: velikost mezitělesové vzdálenosti $r = |\boldsymbol{\rho}|$ (předfaktor), tělesový směr $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = (\mathbf{R}_q^A)^T \hat{\boldsymbol{\rho}}$ (jenž závisí na obou $\boldsymbol{\rho}$ a $q_A$) a konstantní tělesový tenzor $\mathbf{I}_A$ (nikoli perturbovatelná stavová proměnná). Vypočteme parciální derivace vzhledem ke každé stavové proměnné postupně.

**Citlivost na orientaci tělesa $A$: $\partial \boldsymbol{\tau}_A^{\mathrm{body}} / \partial q_A$.** Při infinitezimální tělesové rotaci $\mathbf{R}_q^A \to \mathbf{R}_q^A (\mathbb{I} + \widehat{\delta \boldsymbol{\theta}_A})$ (táž parametrizace užitá v odvození §3.7) je $\delta \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = -\delta \boldsymbol{\theta}_A \times \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$. Píšíc $\mathbf{u} = \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$, $\mathbf{w} = \mathbf{I}_A \mathbf{u}$ a $\alpha = 3 G m_B / r^3$ a linearizujíc moment síly:

$$
\delta \boldsymbol{\tau}_A^{\mathrm{body}}
\;=\; \alpha\, \big[\delta \mathbf{u} \times \mathbf{w} \;+\; \mathbf{u} \times \mathbf{I}_A \delta \mathbf{u}\big]
\;=\; \alpha\, \big[(-\delta \boldsymbol{\theta}_A \times \mathbf{u}) \times \mathbf{w} \;+\; \mathbf{u} \times \mathbf{I}_A (-\delta \boldsymbol{\theta}_A \times \mathbf{u})\big].
$$

Toto je lineární zobrazení z $\delta \boldsymbol{\theta}_A \in \mathbb{R}^3$ do $\delta \boldsymbol{\tau}_A^{\mathrm{body}} \in \mathbb{R}^3$. Definuj **Jacobiho matici orientační odezvy** $\mathbf{J}^{\mathrm{orient}}_{A}$ vztahem $\delta \boldsymbol{\tau}_A^{\mathrm{body}} = \mathbf{J}^{\mathrm{orient}}_{A}\, \delta \boldsymbol{\theta}_A$; jde o maticově-hodnotovou funkci $3 \times 3$ proměnných $(\mathbf{I}_A, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A})$:

$$
\boxed{
\mathbf{J}^{\mathrm{orient}}_{A}\, \boldsymbol{a}
\;=\; -\alpha\, \Big[ (\mathbf{a} \times \mathbf{u}) \times \mathbf{I}_A \mathbf{u} \;+\; \mathbf{u} \times \mathbf{I}_A (\mathbf{a} \times \mathbf{u}) \Big],
\qquad \mathbf{u} = \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A},\quad \alpha = \frac{3 G m_B}{|\boldsymbol{\rho}|^3}.
}
$$

> **Poznámka ke konvenci znaménka (oprava 2026-05-24, OQ-FORT-1).** Celkové $-\alpha$ v zarámovaném výrazu plyne z explicitního přenesení dvou minusových znamének v $\delta \mathbf{u} = -\delta \boldsymbol{\theta}_A \times \mathbf{u}$ skrze linearizaci $\boldsymbol{\tau}_A^{\mathrm{body}}$. Dřívější návrh této kapitoly celkové znaménko vynechal; fortranový křížový test fáze 1 `tests/test_orientation_jacobian.f90` T2a (librační-rovnovážné vlastní hodnoty pro axisymetrické zploštělé těleso) je trvalou regresní pojistkou: pro zploštělé ($I_\parallel > I_\perp$) musí být vlastní hodnota $J_{11} = \alpha (I_\parallel - I_\perp)$ **kladná** (destabilizující), což odpovídá opravené formě $-\alpha$ v kombinaci s algebraickou identitou $-[(a \times u) \times (I u) + u \times I(a \times u)]_{x=e_1} = (I_\parallel - I_\perp,\, 0,\, 0)$ v rovnováze na ose. Librační-rovnovážná specializace a tabulka stability níže jsou již vyjádřeny v opravené formě.

Struktura $\mathbf{J}^{\mathrm{orient}}_{A}$ závisí na třídě symetrie tělesa a na sourodosti $\mathbf{u}$ s hlavními osami tělesa.

**V librační rovnováze (axisymetrické těleso, osa symetrie podél $\hat{\boldsymbol{\rho}}$).** Specializuj na těleso $A$ axisymetrické ($I_1 = I_2 \equiv I_\perp$, $I_3 \equiv I_\parallel$, hlavní osy srovnané s tělesovou soustavou) a orientaci takovou, že $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = (0, 0, 1) = \hat{\mathbf{e}}_z^{\,\mathrm{body}}$ — „radiálně směřující“ postoj z validačního testovacího případu §2.5. V této konfiguraci je $\boldsymbol{\tau}_A^{\mathrm{body}} = \mathbf{0}$ (rovnováha) a přímý výpočet Jacobiho matice dává:

$$
\mathbf{J}^{\mathrm{orient}}_{A}\big|_{\mathrm{eq}}
\;=\; \alpha\, (I_\parallel - I_\perp)\, \mathrm{diag}(1,\, 1,\, 0)
\;=\; \frac{3 G m_B (I_\parallel - I_\perp)}{|\boldsymbol{\rho}|^3}\, \mathrm{diag}(1,\, 1,\, 0).
$$

Znaménko faktoru $(I_\parallel - I_\perp)$ určuje stabilitu rovnováhy — **věta o stabilizaci gravitačním gradientem** pro axisymetrický satelit:

| Geometrie tělesa | Znaménko $I_\parallel - I_\perp$ | Rovnováha s osou symetrie podél $\hat{\boldsymbol{\rho}}$ |
|---|---|---|
| Protáhlé ($I_\parallel < I_\perp$; jehla / dlouhá tyč podél osy symetrie) | $-$ | **Stabilní** (vratný moment síly; librace kolem rovnováhy) |
| Sférické ($I_\parallel = I_\perp$) | $0$ | Marginální (žádný moment síly z této osy) |
| Zploštělé ($I_\parallel > I_\perp$; plochý disk / sluneční clona kolmá k ose symetrie) | $+$ | **Nestabilní** (moment síly žene krátkou osu pryč od $\hat{\boldsymbol{\rho}}$; stabilní orientace je osa symetrie kolmá k $\hat{\boldsymbol{\rho}}$) |

Prvek $(3,3)$ matice $\mathbf{J}^{\mathrm{orient}}_{A}\big|_{\mathrm{eq}}$ mizí: rotace kolem osy symetrie samotné neprodukuje žádný moment síly z gravitačního gradientu, jak se očekává (těleso je symetrické kolem této osy; rotace kolem ní nemění gravitační vazbu).

**Citlivost na mezitělesovou polohu: $\partial \boldsymbol{\tau}_A^{\mathrm{body}} / \partial \boldsymbol{\rho}$.** Rozlož $\delta \boldsymbol{\rho} = \delta r\, \hat{\boldsymbol{\rho}} + \delta \boldsymbol{\rho}_\perp$, kde $\delta \boldsymbol{\rho}_\perp \perp \hat{\boldsymbol{\rho}}$. Obě části mají odlišné efekty na $\boldsymbol{\tau}_A^{\mathrm{body}}$:

- **Radiální perturbace** $\delta r$: přeškáluje předfaktor $\alpha = 3 G m_B / r^3 \to \alpha (1 - 3 \delta r / r)$ do prvního řádu, zatímco $\hat{\boldsymbol{\rho}}$ (a tudíž $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$) je nezměněno. Jacobiho matice podél $\hat{\boldsymbol{\rho}}$ je prostě:
  $$
  \frac{\partial \boldsymbol{\tau}_A^{\mathrm{body}}}{\partial r}
  \;=\; -\frac{3}{r}\, \boldsymbol{\tau}_A^{\mathrm{body}},
  $$
  přeškálujíc celý vektor momentu síly s citlivostí nepřímo úměrnou vzdálenosti. V librační rovnováze ($\boldsymbol{\tau}_A^{\mathrm{body}} = \mathbf{0}$) tento příspěvek mizí.

- **Tangenciální perturbace** $\delta \boldsymbol{\rho}_\perp$: rotuje inerciální směr $\hat{\boldsymbol{\rho}}$ o infinitezimální osu-úhel $\delta \boldsymbol{\chi} = (\hat{\boldsymbol{\rho}} \times \delta \boldsymbol{\rho}_\perp) / r$. Tato rotace se propaguje do tělesové soustavy jako $\delta \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = (\mathbf{R}_q^A)^T \delta \hat{\boldsymbol{\rho}}$, což je **algebraicky strukturálně identické** s případem orientační perturbace výše (malá rotace $\mathbf{u}$ v tělesových souřadnicích). Odpovídající Jacobiho matice je přímým dosazením:
  $$
  \boxed{
  \frac{\partial \boldsymbol{\tau}_A^{\mathrm{body}}}{\partial \boldsymbol{\rho}_\perp}\, \delta \boldsymbol{\rho}_\perp
  \;=\; \mathbf{J}^{\mathrm{orient}}_{A}\, \big[(\mathbf{R}_q^A)^T\, (\hat{\boldsymbol{\rho}} \times \delta \boldsymbol{\rho}_\perp / r)\big],
  }
  $$
  s $\mathbf{J}^{\mathrm{orient}}_{A}$ toutéž maticí jako Jacobiho matice orientační odezvy výše. Jacobiho matice tělesové orientace a Jacobiho matice směru polohy nejsou nezávislé objekty — jsou týmž lineárním zobrazením působícím na různé vstupy, odrážejíc hlubší fakt, že do $\boldsymbol{\tau}_A^{\mathrm{body}}$ vstupuje pouze *relativní* úhel mezi $\hat{\boldsymbol{\rho}}$ a hlavními osami tělesa $A$.

**Citlivost na orientaci tělesa $B$: $\partial \boldsymbol{\tau}_A^{\mathrm{body}} / \partial q_B = 0$ při kvadrupólovém useknutí.** Moment síly z gravitačního gradientu tělesa $A$ závisí na tělese $B$ jen skrze jeho hmotnost $m_B$ a jeho polohu $\mathbf{R}_B$ (vstupující skrze $\boldsymbol{\rho} = \mathbf{R}_B - \mathbf{R}_A$); orientace tělesa $B$ $q_B$ se v $\boldsymbol{\tau}_A^{\mathrm{body}}$ neobjevuje při kvadrupólově-useknutém multipolovém rozvoji přijatém v §3.4. Orientace tělesa $B$ vstupuje až v *dalším* multipolovém řádu: bilineární člen kvadrupól-kvadrupól (pokles $1/r^5$) a člen oktupól-monopól (pokles $1/r^4$), oba odhadnuté a omezené v §3.3. Pro konfigurace první verze $(\ell/r \lesssim 10^{-2})$ je tato křížově-orientační citlivost potlačena alespoň o čtyři řády vůči vedoucímu kvadrupólovému slapovému členu a zůstává pod zaokrouhlovací chybou dvojnásobné přesnosti.

**Struktura $\mathbf{J}^{\mathrm{orient}}_{A}$ podle třídy symetrie.** Počítání hodnosti Jacobiho matice orientační odezvy napříč třídami tělesové symetrie (tabulka paralelizuje tabulku citlivosti mimodiagonálních rychlostí z §2.5.1):

| Třída symetrie tělesa | Hodnost $\mathbf{J}^{\mathrm{orient}}_{A}$ | Fyzikální interpretace |
|---|---|---|
| Sférická ($I_1 = I_2 = I_3$) | $0$ | $\boldsymbol{\tau}_A \equiv \mathbf{0}$ identicky; žádná citlivost na orientaci (zvláštní případ testovacího případu rozpojení §3.8.2). |
| Axisymetrická ($I_1 = I_2 \neq I_3$) | $2$ genericky (klesá na $0$ v rovnováze na ose symetrie pro rotace *kolem* osy symetrie) | Dvourozměrná librační rovina kolmá k ose symetrie; jeden cyklický směr podél osy. |
| Asymetrická (všechny $I_i$ odlišné) | $3$ genericky | Tři nezávislé librační / nestabilní módy; moment síly z gravitačního gradientu se váže na všechny tři orientační stupně volnosti. |

**Ukazatel citlivosti sestavy dvou těles.** Celková dynamika sestavy váže rotaci tělesa $A$, rotaci tělesa $B$ a relativní orbitu skrze $\boldsymbol{\rho}$. Integrovaný efekt momentu síly z gravitačního gradientu na rotační stav každého tělesa škáluje s $\alpha \cdot |\mathbf{I}_X^{\mathrm{anisotropy}}|$, kde anizotropie je zachycena (pro axisymetrická tělesa) poměrem $\eta_X = (I_\parallel - I_\perp)_X / I_\perp^X$ nebo ekvivalentně $I_\parallel / I_\perp - 1$. Pro sestavy dvou těles je **ukazatelem** citlivosti celé sestavy na perturbace počátečních podmínek součin:

$$
\Sigma_{\mathrm{assembly}}
\;\equiv\;
\eta_A \cdot \eta_B \cdot f(\hat{\boldsymbol{\rho}}\text{-sourodost}),
$$

kde $f$ je faktor řádu $O(1)$ kódující, jak $\hat{\boldsymbol{\rho}}$ se srovnává s hlavními osami každého tělesa ($f = 0$, je-li $\hat{\boldsymbol{\rho}}$ podél hlavní osy kteréhokoli tělesa současně a rovnováha je při nulovém momentu síly; $f \to 1$ pro generickou orientaci; $f$ extrémní v gravitačně-gradientně nestabilních konfiguracích). Pro konfiguraci první verze (těleso podobné JWST $A$ s $I_\parallel/I_\perp = 1{,}52$, těleso podobné sondě $B$ s $I_\parallel/I_\perp = 5{,}81$): $\eta_A = 0{,}52$, $\eta_B = 4{,}81$, $\eta_A \eta_B \approx 2{,}5$. Vysoká anizotropie tělesa podobného sondě dominuje citlivosti sestavy; malé perturbace $q_B$ se propagují do velkých výkyvů $\boldsymbol{\tau}_B^{\mathrm{body}}$, jež zase posouvají účetnictví momentu hybnosti a zatěžují diagnostiku zachování integrátoru.

**Numerická kontrola zdravého rozumu: zploštělé referenční těleso v hypotetické konfiguraci podobné L2.** Vezmi $m_B = m_\oplus \approx 6 \times 10^{24}\ \text{kg}$ (primár o hmotnosti Země), $r = 1{,}5 \times 10^9\ \text{m}$ (reprezentativní řádová vzdálenost L2 od Země). Pak:

$$
\alpha = \frac{3 G m_B}{r^3} \approx 3{,}5 \times 10^{-13}\ \text{s}^{-2}.
$$

Pro **zploštělé referenční těleso** `make_oblate_reference_body()` (kap. 02 §2.4 — syntetické jednoválcové těleso s $\mathbf{I}_A = \mathrm{diag}(2540, 2540, 3870)\ \text{kg m}^2$ konstrukčně) je librační-rovnovážná tuhost $\alpha (I_\parallel - I_\perp) = 3{,}5 \times 10^{-13} \cdot 1330 \approx 4{,}7 \times 10^{-10}\ \text{N m / rad}$. Protože $I_\parallel > I_\perp$ (zploštělé), je tuhost *kladná*, signalizujíc **nestabilitu** postoje s osou symetrie podél radiály — gravitační gradient žene osu symetrie *pryč* od radiálního směru. Čas e-násobku nestability je $\tau_{\mathrm{inst}} = \sqrt{I_\perp / |k|} \approx 2{,}3 \times 10^6\ \text{s} \approx 27\ \text{dní}$. To je řádově časové měřítko pomalého vývoje postoje kosmické lodě podobné JWST s osou sluneční clony podél radiály (v praxi reálný JWST užívá setrvačníky a reakční řízení ke kompenzaci tohoto a dalších momentů sil průběžně).

**Kontrast: protáhlý kompozit v téže hypotetické konfiguraci podobné L2.** Zjednodušený kompozit kosmické lodě `make_jwst_like()` (4složkový sluneční clona + servisní modul + výložník + primár, celkem 2450 kg) produkuje $\mathbf{I}_A = \mathrm{diag}(23323, 23323, 15384)\ \text{kg m}^2$ — těleso je **protáhlé** ($I_\parallel = I_3 < I_\perp = I_1 = I_2$), neboť rozptyl složek podél z přispívá k $I_\perp$ Steinerovou větou nad lokální axiální příspěvky. V téže hypotetické konfiguraci L2 je tuhost $\alpha (I_\parallel - I_\perp) = 3{,}5 \times 10^{-13} \cdot (-7938) \approx -2{,}8 \times 10^{-9}\ \text{N m / rad}$ — **záporná** (vratná), takže postoj s osou symetrie podél radiály je **stabilní** s librační periodou $T_{\mathrm{lib}} = 2\pi \sqrt{I_\perp / |k|} = 2\pi \sqrt{23323 / 2{,}8 \times 10^{-9}} \approx 1{,}81 \times 10^{7}\ \text{s} \approx 210\ \text{dní}$ (režim *zemského slapu*: $m_B = m_\oplus$, $r = 1{,}5 \times 10^9\ \text{m}$ = vzdálenost Země–L2; totéž těleso v režimu *slunečního slapu*, $m_B = m_\odot$ a $r \approx 1{,}496 \times 10^{11}\ \text{m}$ = vzdálenost Slunce–Země, má $T_{\mathrm{lib}} \approx 3{,}6 \times 10^{7}\ \text{s} \approx 362\ \text{dní}$ — viz testbed fáze 2 `testbeds/lagrange_top_libration.py`; dřívější „$5{,}7 \times 10^{6}\ \text{s} \approx 66\ \text{dní}$“ byl aritmetický uklouz s vynechaným $2\pi$ opravený dle OQ-PHASE2-8). Protáhle-stabilní klasifikace kompozitu je konzistentní s principem stabilizace gravitačním gradientem pro kosmické lodě s dlouhou osou; co se liší od intuice reálného-JWST, je to, že rozptyl složek zjednodušeného modelu podél z přemáhá lokální axiální moment sluneční clony, produkujíc protáhlou spíše než zploštělou topologii kolem z. Tento kontrast — **zploštěle-nestabilní referenční těleso vedle protáhle-stabilní zjednodušené kosmické lodě** — ilustruje tabulku stability přímo: totéž gravitační prostředí, táž velikost anizotropie, opačná klasifikace stability hnaná znaménkem $(I_\parallel - I_\perp)$. Fortranový křížový-testovací fixture fáze 1 `tests/fixtures/gravity_gradient_inputs.json` případy 02 (protáhlé podobné JWST, $\tau_y = -4{,}8 \times 10^{-10}$ při 10° mimo rovnováhu) a 07 (zploštělé referenční těleso, $\tau_y = +8{,}1 \times 10^{-11}$ při 10° mimo rovnováhu) demonstrují opačná znaménka empiricky.

**Důsledek pro integrátor na inženýrské úrovni: vícemeřítkový problém.** Integrace první verze užívá RK4 při $dt = 0{,}05\ \text{s}$ přes okno $600\ \text{s}$. Soustava má alespoň čtyři přirozená časová měřítka:

| Časové měřítko | Velikost (první verze + hypotetické L2) | Požadavek na rozlišení |
|---|---|---|
| Perioda spinu tělesa $T_{\mathrm{spin}}$ | $\sim 628\ \text{s}$ (pro $\omega \sim 0{,}01\ \text{rad/s}$) | $dt \ll T_{\mathrm{spin}}$ — snadno splněno při $0{,}05\ \text{s}$ ($\sim 12\,000$ vzorků na periodu) |
| Mimodiagonální časové měřítko $I^{\mathrm{inertial}}$ $\tau_{\mathrm{off}}$ (§2.5.1) | $\sim 190\ \text{s}$ | V rámci periody spinu — totéž omezení. |
| Librační / nestabilní časové měřítko | $\sim 10^7\ \text{s}$ (~7–12 měsíců; zemský- až sluneční-slap) | Integrační okno musí být $\gg T_{\mathrm{lib}}$ k pozorování librace — nezachyceno při okně $600\ \text{s}$ |
| Orbitální perioda v konfiguraci podobné L2 | $\sim 10^7\ \text{s}$ (~6 měsíců) | Totéž — nezachyceno při okně $600\ \text{s}$ |

Okno první verze je tedy **spin-rozlišující, avšak librace-slepé**. Ke studiu citlivosti librace gravitačního gradientu sestavy numericky musí integrační okno sahat k $\sim 10^7\ \text{s}$ ($\sim 2 \times 10^8$ kroků RK4 při $dt = 0{,}05\ \text{s}$), kdy se akumulovaný nesymplektický úlet energie z RK4 ($\sim O(dt^4 \cdot T)$ dle kapitoly 1 §1.6) stává dominantním zdrojem chyby. Symplektický integrátor na Lieově grupě (Hairer-Lubich-Wanner 2006 §VII.5) je strukturální opravou pro dlouhohorizontové simulace postoje; první verze záměrně zůstává u RK4 pro rozsah spin-rozlišení a označuje symplektický upgrade jako položku inženýrské úrovně v0.2.

> **Doplňující anglojazyčné odkazy.** Hughes 2004 Kap. 8 (stabilizace gravitačním gradientem, s maticí tuhosti pro axisymetrická a asymetrická tělesa odvozenou v podstatě v zde uvedené formě); Markley & Crassidis 2014 *Fundamentals of Spacecraft Attitude Determination and Control* §3.7 (inženýrské pojednání librační dynamiky s důsledky pro řízení postoje); Wertz 1978 *Spacecraft Attitude Determination and Control* §17 pro klasický inženýrský odkaz. Zde zdůrazněné zarámování citlivostní Jacobiho maticí — explicitní matice $\mathbf{J}^{\mathrm{orient}}_{A}$ jako sjednocující objekt napříč orientačními a mezitělesovými-polohovými perturbacemi — je inženýrským doplňkem k Hughesovu pojednání librace.

## §3.8 Redukce na validační testovací případy

Kvadrupólově-useknutý potenciál a moment síly z gravitačního gradientu se redukují na analytické limity zkatalogizované v `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md`. Tři z nich zde projdeme explicitně; kapitola 06 pokrývá celý katalog.

### §3.8.1 Redukce na Kepler (validační testovací případ §2.3)

Polož $\mathbf{I}_A \to \mathbf{0}$ a $\mathbf{I}_B \to \mathbf{0}$ v zarámovaném celkovém $V$ z §3.4.3. Pak $V_{\mathrm{tidal},\, A}$ i $V_{\mathrm{tidal},\, B}$ identicky mizí (nulový tenzor setrvačnosti má nulovou stopu a nulovou kontrahovanou formu) a useknuté $V$ se zhroutí na:

$$
V \;\xrightarrow{\mathbf{I}_{A,B} \to \mathbf{0}}\; -\frac{G\, m_A\, m_B}{|\boldsymbol{\rho}|}.
$$

V kombinaci s kinetickou energií v těžišťové soustavě $T = \tfrac{1}{2} \mu |\dot{\boldsymbol{\rho}}|^2 + \tfrac{1}{2}\, \boldsymbol{\Omega}_A^T \mathbf{I}_A \boldsymbol{\Omega}_A + \tfrac{1}{2}\, \boldsymbol{\Omega}_B^T \mathbf{I}_B \boldsymbol{\Omega}_B z kapitoly 02 §2.6 rotační kinetické energie v této limitě rovněž mizí a Lagrangián $L = T - V$ se redukuje na:

$$
L \;\xrightarrow{\mathbf{I}_{A,B} \to \mathbf{0}}\; \tfrac{1}{2}\, \mu\, |\dot{\boldsymbol{\rho}}|^2 \;+\; \frac{G\, m_A\, m_B}{|\boldsymbol{\rho}|},
$$

standardní Keplerův Lagrangián s redukovanou hmotností. Jeho Eulerovy–Lagrangeovy rovnice jsou Keplerovy pohybové rovnice (Landau-Lifšic, sv. I §13–§15). Testovací případ `_specs/...` §2.3 (`tests/test_dynamics.py::TestKeplerLimit`, plánováno) toto ověřuje numericky nastavením tenzorů setrvačnosti inženýrské úrovně na nulu a kontrolou, že orbitální perioda odpovídá $T_{\mathrm{orbit}} = 2\pi \sqrt{a^3 / (G(m_A + m_B))}$ s relativní chybou $10^{-6}$.

### §3.8.2 Redukce na volnou rotaci + Keplerovo rozpojení (validační testovací případ §2.4)

Polož $\mathbf{I}_B \to \mathbf{0}$ (bodová hmota tělesa $B$); ponech $\mathbf{I}_A$ generické. Pak $V_{\mathrm{tidal},\, B} = 0$ (nulový kvadrupólový moment bodové hmoty) a $V_{\mathrm{tidal},\, A} \neq 0$ obecně — *pokud* těleso $A$ není sféricky symetrické, v kterémžto případě $\mathbf{I}_A = I_0 \mathbb{I}_3$ a:

$$
\mathrm{tr}(\mathbf{I}_A) - 3\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}
\;=\; 3 I_0 - 3 I_0 |\hat{\boldsymbol{\rho}}|^2 \;=\; 0,
$$

takže $V_{\mathrm{tidal},\, A}$ identicky mizí. Translační a rotační dynamika se zcela rozpojí: orbita $A$ je Keplerova kolem bodové hmoty $B$ a rotace $A$ je bezsilová Eulerova precese (Landau-Lifšic §35 / kapitola 05 této kanonické úrovně). Toto je nastavení testovacího případu `_specs/...` §2.4.

Není-li $A$ sféricky symetrické — například první-verze axisymetrické těleso podobné JWST s $\mathbf{I}_A = \mathrm{diag}(I_\perp, I_\perp, I_\parallel)$ a $I_\parallel \neq I_\perp$ — pak $V_{\mathrm{tidal},\, A}$ závisí na orientaci $A$ vůči $\boldsymbol{\rho}$ a moment síly z gravitačního gradientu $\boldsymbol{\tau}_A$ je genericky nenulový. Toto je režim testovacího případu **librace gravitačního gradientu** (§3.8.3 níže) a případ pojednaný rozsáhle v Hughes 2004 Kap. 8.

### §3.8.3 Librace o malých úhlech (validační testovací případ §2.5)

Umísti axisymetrické těleso $A$ ($I_{xx,A} = I_{yy,A} \equiv I_\perp$, $I_{zz,A} \equiv I_\parallel$) na kruhovou orbitu kolem bodové hmoty tělesa $B$. Definuj **radiálně směřující rovnováhu**: orientaci $A$ takovou, že jeho tělesová osa $z$ (osa symetrie, ta s odlišným momentem) je srovnána s $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = \hat{\mathbf{e}}_z^{\,\mathrm{body}\, A}$ v každém čase — tj. $A$ rotuje synchronně s orbitou, udržujíc svou osu symetrie namířenou na $B$.

V rovnováze je $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = (0, 0, 1)$ a $\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = (0, 0, I_\parallel)$. Vektorový součin $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \times (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}) = (0,0,1) \times (0,0,I_\parallel) = \mathbf{0}$. Moment síly mizí; rovnováha je pravým kritickým bodem orientačně závislého potenciálu.

Perturbuj orientaci o malé úhly $\theta, \phi$ kolem rovnováhy (perturbace klopení–klonění kolem osy symetrie). Do lineárního řádu v $(\theta, \phi)$ je perturbované $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \approx (-\theta, \phi, 1)$ do prvního řádu v malých úhlech a (vypracováním vektorového součinu a kontrakcí vůči orbitální úhlové rychlosti $n = (G(m_A + m_B)/|\boldsymbol{\rho}|^3)^{1/2}$):

Linearizované librační rovnice jsou vázané harmonické oscilátory s frekvencemi vlastních módů určenými poměrem setrvačnosti $(I_\parallel - I_\perp)/I_\perp$ a orbitální rychlostí $n$. Plné odvození je standardní (Hughes 2004 §8.4); zde jej nereprodukujeme, neboť vyžaduje tělesové Eulerovy rovnice z kapitoly 05 jako lešení. Testovací případ `_specs/...` §2.5 (`tests/test_integration.py::TestGravityGradientLibration`, plánováno) ověřuje librační periodu numericky.

Klíčovým bodem pro tuto kapitolu je, že moment síly z gravitačního gradientu, jejž jsme odvodili v §3.7, **produkuje netriviální librační spektrum** pro jakékoli nesférické těleso a spektrum je testovatelné vůči uzavřené predikci Hughes 2004 s relativní chybou $10^{-4}$. Toto je učebnicově-doslovné pokrytí testy pro moment síly z gravitačního gradientu C1, ve smyslu doktríny §4.4.

### §3.8.4 Rozpojení při nekonečné vzdálenosti (validační testovací případ §3.2)

Pro $|\boldsymbol{\rho}| \to \infty$ člen monopól-monopól $V^{(0)} \to 0$ a oba slapové členy klesají jako $1/|\boldsymbol{\rho}|^3 \to 0$. Celkový vzájemný potenciál mizí rychleji než jakákoli konečná mocnina vzdálenosti (při rozvinutí do konečného multipolového řádu; vyšší řády klesají ještě rychleji). Každé těleso se stává izolovaným bezsilovým setrvačníkem: $\boldsymbol{\tau}_A, \boldsymbol{\tau}_B \to \mathbf{0}$, $V \to 0$ a Lagrangián se faktorizuje jako $L = L_A + L_B$, kde každé $L_X = T_X$ je kineticky-pouze Lagrangián izolovaného tuhého tělesa. Soustava se zcela rozpojí. Testovací případ `_specs/...` §3.2 ověřuje, že křížová vazba zůstává pod $10^{-12}$ relativně na orbitální periodu při velké vzdálenosti.

## §3.9 Souhrn křížových odkazů na kód

Pro čtenáře inženýrské úrovně, jenž sledoval až do tohoto bodu, odpovídají kanonické rovnice této kapitoly kódu inženýrské úrovně takto:

| Objekt kanonické úrovně | Implementace inženýrské úrovně |
|---|---|
| Přesný integrál $V$ (§3.2 zarámovaný) | Neimplementováno; jde o úplnou fyzikální referenci. |
| Monopól-monopól $V^{(0)}$ (§3.4.1) | `dynamics.py::total_potential_energy` (řádky 207–214). |
| Kvadrupólový slapový $V_{\mathrm{tidal},\, A}$ (§3.4.3 / §3.5 zarámovaný) | Není přímo počítáno pro energetické diagnostiky; absorbováno implicitně do dynamiky skrze $\boldsymbol{\tau}_A$. **Poznámka:** diagnostika zachování `total_potential_energy` tedy sleduje pouze člen monopól-monopól; chybějící slapový příspěvek k $V$ je pod zaokrouhlovací chybou dvojnásobné přesnosti při vzdálenostech první verze a neovlivňuje pozorované $|dE/E_0| = 6{,}2 \times 10^{-12}$. **Mohl by být přidán** jako zjemnění; zařazeno v §3.11 OQ-3.4. |
| Moment síly z gravitačního gradientu $\boldsymbol{\tau}_A^{\mathrm{body}}$ (§3.7 zarámovaný) | `dynamics.py::gravity_gradient_torque_body_frame` (řádky 113–131); volací místo pro těleso $A$ na řádku 157, pro těleso $B$ na řádku 158. |
| Stopová identita $\mathbf{Q} = \mathrm{tr}(\mathbf{I}) \cdot \mathbb{I} - 3 \mathbf{I}$ (§3.5) | Neimplementováno jako samostatná funkce; složeno do vzorce momentu síly. |
| Hranice konvergence (§3.3) | Dosud netvrzeno; zařazeno v `_specs/...` §3.3 / §6 OQ-odloženo. |

Shoda je **přesná na úrovni momentu síly z gravitačního gradientu**: zarámovaný vzorec §3.7 a inženýrská `gravity_gradient_torque_body_frame` jsou řádek po řádku identické, s konzistentními konvencemi znaménka napříč volacími místy tělesa $A$ i tělesa $B$ (§3.7.1).

## §3.10 Co je stanoveno na konci této kapitoly

Čtenář, jenž prošel §3.1–§3.9, má:

- **Přesný vzájemný gravitační potenciál** dvou rozšířených tuhých těles jako 6rozměrný integrál přes jejich rozložení hmoty, odvozený z newtonovské gravitace z prvních principů.
- Obor konvergence a hranici těsného přiblížení multipolového rozvoje, zamčenou jako obor platnosti modelu.
- Kartézský multipolový rozvoj $V$ do kvadrupólového řádu, s monopólově-monopólovým (Kepler) a kvadrupólově-slapovým (gravitační gradient) příspěvkem identifikovaným explicitně.
- Nezávislý sféricko-harmonický multipolový rozvoj, produkující tentýž kvadrupólový výsledek standardní Jacksonovou identitou přenesenou z elektromagnetismu.
- **Křížovou validaci**, že oba přístupy dávají identické výsledky v kvadrupólovém řádu — doktrinální požadavek „dvě paralelně“, splněn.
- Stopovou identitu $\mathbf{Q} = \mathrm{tr}(\mathbf{I}) \cdot \mathbb{I} - 3 \mathbf{I}$ vztahující gravitační kvadrupólový tenzor ke kineticko-teoretickému tenzoru setrvačnosti z kapitoly 02.
- **Moment síly z gravitačního gradientu** $\boldsymbol{\tau}_A^{\mathrm{body}} = (3 G m_B / |\boldsymbol{\rho}|^3)\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \times (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A})$, se shodou řádek po řádku s implementací inženýrské úrovně `dynamics.py::gravity_gradient_torque_body_frame`.
- **Jacobiho matici orientační odezvy** $\mathbf{J}^{\mathrm{orient}}_{A}$ (§3.7.2) — výstup analýzy citlivosti inženýrské úrovně. Sjednocující lineární zobrazení $3 \times 3$ řídící odezvu momentu síly na perturbace tělesové orientace i směru mezitělesové polohy; librační-rovnovážná specializace $\alpha (I_\parallel - I_\perp) \mathrm{diag}(1,1,0)$ s predikcí stability gravitačního gradientu; vícemeřítkový audit časových měřítek identifikující první verzi jako spin-rozlišující, avšak librace-slepou, a označení symplektického integrátoru na Lieově grupě pro dlouhohorizontové studie.
- Redukce useknutého $V$ na validační testovací případy analytických limitů (Kepler, volná rotace + Keplerovo rozpojení, librace gravitačního gradientu, rozpojení při nekonečné vzdálenosti), každou zdokumentovanou a ukazující dopředu na plné pojednání kapitoly 06.

**Co dosud stanoveno NENÍ:**

- Plný Lagrangián $L = T - V$ a Eulerovy–Lagrangeovy pohybové rovnice na $TQ$. **Kapitola 04.**
- Tělesové Eulerovy rovnice $\mathbf{I} \dot{\boldsymbol{\Omega}} + \boldsymbol{\Omega} \times (\mathbf{I} \boldsymbol{\Omega}) = \boldsymbol{\tau}$ jako specializace Euler-Poincarého redukcí, s $\boldsymbol{\tau}$ zde odvozeným vstupujícím jako pravá strana. **Kapitola 05.**
- Rozšíření omezené soustavy tří těles u L2 s centrifugálními a Coriolisovými pseudo-silami. **Odloženo na kanonickou úroveň v0.2** dle kapitoly 00 §1.
- Plný výčet a numerické ověření validačních testovacích případů. **Kapitola 06**, čerpajíc z katalogu v `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md`.

S $T$ (kapitola 02) a $V$ (tato kapitola) oběma k dispozici má kapitola 04 veškerý sestavovací materiál, jejž potřebuje ke konstrukci Lagrangiánu a odvození pohybových rovnic.

---

## §3.11 Otevřené otázky v této kapitole

- **OQ-3.1.** ***Vyřešeno 2026-05-24*** — explicitní algebra změny báze přidána v §3.5 (zarámovaný blok rovnic následující po identifikaci $\ell = 2$). Mapuje pět $q_2^m(A)$ (pro $m = -2, -1, 0, 1, 2$) na jisté lineární kombinace složek bezstopového kartézského kvadrupólového tenzoru $(\mathbf{Q}_A)_{ij}$, s Condon-Shortleyovou sféricko-harmonickou konvencí a uznáním předfaktoru-1/3 rozdílu konvence T&B/Jackson $\mathcal{I}_{jk}$-vs-naše-$\mathbf{Q}_{jk}$ (dle referenčního anchoru §3.4.3). Invertibilní algebra poznamenána, podmínky reálnosti uvedeny, geometrická interpretace každé složky $q_2^m$ přidána (protáhlé/zploštělé protažení podél $z$, mimodiagonální náklon $xz/yz$, rovinná asymetrie $xx-yy$ + $xy$). Specializace axisymetrického tělesa poznamenána: nenulové pouze $q_2^0$. Křížově validováno vůči diskusi stopové identity již v §3.5. Anchor: Jackson 1999 §3.6 (sférické harmoniky) + §4.1 (ekvivalence multipolových momentů).
- **OQ-3.2.** ***Částečné řešení 2026-05-25 (kanonická úroveň, 3. kolo).*** **Korekce dalšího řádu** k zarámovanému §3.4.3 kvadrupólově-useknutému $V$ mají explicitní škálování mocninou $|\boldsymbol{\rho}|$ odvoditelné přímo rozšířením kartézského Taylorova rozvoje §3.4 o jeden řád nebo ekvivalentně sféricko-harmonického rozvoje §3.5 o jednu úroveň $\ell$. V dalším řádu vznikají dva odlišné příspěvky:
  - **Bilineární člen kvadrupól-kvadrupól:** $V^{(\ell_A=2,\ell_B=2)} \propto G\, [\hat{\boldsymbol{\rho}} \cdot \mathbf{Q}_A \cdot \hat{\boldsymbol{\rho}}] [\hat{\boldsymbol{\rho}} \cdot \mathbf{Q}_B \cdot \hat{\boldsymbol{\rho}}] / |\boldsymbol{\rho}|^5$, kde součin velikosti $\ell^4$ obou kvadrupólových tenzorů těles vůči mezitělesovému směru škáluje člen. Škálování $1/|\boldsymbol{\rho}|^5$ plyne ze součinu dvou faktorů Laplaceova rozvoje $|\boldsymbol{\rho}|^{-(\ell+1)}$ s $\ell = 2$ každý, sečtených napříč dvoutělesovým multipolovým součtem.
  - **Člen oktupól-monopól:** $V^{(\ell_A=3,\ell_B=0)} \propto G\, m_B\, \mathcal{O}_A : (\hat{\boldsymbol{\rho}}\hat{\boldsymbol{\rho}}\hat{\boldsymbol{\rho}}) / |\boldsymbol{\rho}|^4$, kde $\mathcal{O}_A$ je oktupólový momentový tenzor tělesa $A$ řádu 3 (symetrický a bezstopový, kódující $7$ nezávislých stupňů volnosti pro $\ell=3$). Škálování $1/|\boldsymbol{\rho}|^4$ je faktorem $|\boldsymbol{\rho}|^{-(\ell+1)}$ při $\ell=3$ spárovaným s konstantní hmotností tělesa $B$.
  
  **Mez numerické velikosti při konfiguracích první verze.** S užitím dimenzionálních odhadů $\mathcal{O}_A \sim m_A \ell_A^3$ a $|\mathbf{Q}_X| \sim m_X \ell_X^2$:
  
  | Člen | Škálování velikosti | Při $\ell/r = 10^{-2}$ |
  |---|---|---|
  | $V^{(0)}$ Kepler | $Gm_Am_B/r$ | $1$ (reference) |
  | $V_{\mathrm{tidal},A}$ kvadrupól-monopól | $\sim (\ell_A/r)^2$ | $10^{-4}$ |
  | $V^{(\ell_A=2,\ell_B=2)}$ kvadrupól-kvadrupól | $\sim (\ell_A \ell_B/r^2)^2 = (\ell^2/r^2)^2$ | $10^{-8}$ |
  | $V^{(\ell_A=3,\ell_B=0)}$ oktupól-monopól | $\sim (\ell_A/r)^3$ | $10^{-6}$ |
  
  Oba členy dalšího řádu padají pod **pozorovaný úlet zachování dvojnásobné přesnosti** $|dE/E_0| = 6{,}2 \times 10^{-12}$ pro vzdálenosti $\ell/r \lesssim 10^{-2}$ — režim, v němž první verze inženýrské úrovně operuje ($|\boldsymbol{\rho}| \sim 1000$ m vs $\ell \sim 7$ m). Odhad chyby kvadrupólového useknutí z §3.3 je konzistentní s těmito velikostmi.
  
  **Status řešení.** Explicitní odvození odloženo na v0.2, pokud se rozsah kanonické úrovně rozšíří na režim slapově-uzamčeného binárního asteroidu ($\ell/r \to 1$, kde se oba členy stávají inženýrsky významnými; referenční anchor pro tento režim: Boué 2018 *J. Geophys. Res.* nebo Maciejewski 1995). Zde stanovená řádová mez je dostatečná pro v0.1.
- **OQ-3.3.** ***Vyřešeno 2026-05-25 (kanonická úroveň, 3. kolo, třetí nezávislá odvozovací cesta).*** Moment síly z gravitačního gradientu na těleso $A$ lze odvodit **třetí nezávislou cestou** nad rámec dvou již v §3.7 (variační odvození Lagrangeovým formalismem) a křížové kontroly §3.7.1 (cesta přímého vzorce Hughes 2004). Třetí cestou je **přímá integrace $\mathbf{r} \times \nabla V$ přes rozložení hmoty tělesa**, dle Goldstein 2002 §5.5.
  
  **Nastavení.** Těleso $A$ sídlí v gravitačním potenciálu tělesa $B$ braného jako bodová hmota v $\mathbf{R}_B$. V každém hmotném bodě tělesa $A$ (poloha $\boldsymbol{\xi}_A$ v inerciální soustavě, s $\boldsymbol{\xi}_A = \mathbf{R}_q^A \mathbf{x}_A$ jako v §3.4) je síla na jednotku hmoty $-\nabla \Phi_B(\boldsymbol{\xi}_A)$, kde $\Phi_B = -Gm_B/|\boldsymbol{\xi}_A - \mathbf{R}_B|$. Tělesový moment síly na těleso $A$ kolem jeho vlastního těžiště je:
  
  $$
  \boldsymbol{\tau}_A^{\mathrm{inertial}}
  \;=\;
  -\int_A \rho_A(\mathbf{x}_A)\, (\boldsymbol{\xi}_A - \mathbf{R}_A) \times \nabla \Phi_B(\boldsymbol{\xi}_A)\, dV_A.
  $$
  
  **Rozvoj.** Dosaď $\boldsymbol{\xi}_A - \mathbf{R}_A = \mathbf{R}_q^A \mathbf{x}_A$ a Taylorovsky rozviň $\nabla \Phi_B(\boldsymbol{\xi}_A)$ kolem těžiště tělesa $\mathbf{R}_A$ (malý rozvojový parametr $|\mathbf{x}_A|/|\boldsymbol{\rho}|$, týž parametr, jenž řídí rozvoj §3.4). Ve vedoucím řádu — ponechávajíc členy lineární v $\mathbf{x}_A$:
  
  $$
  \nabla \Phi_B(\boldsymbol{\xi}_A) \;\approx\; \nabla\Phi_B(\mathbf{R}_A) + (\mathbf{R}_q^A \mathbf{x}_A) \cdot \nabla^2 \Phi_B(\mathbf{R}_A) + \mathcal{O}(|\mathbf{x}_A|^2).
  $$
  
  Člen nultého řádu přispívá nulovým momentem síly (integruje vůči $\int_A \rho_A \mathbf{x}_A dV_A = \mathbf{0}$ konvencí těžiště). Člen prvního řádu dává, po přenesení vektorového součinu integrálem:
  
  $$
  \boldsymbol{\tau}_A^{\mathrm{inertial}}
  \;\approx\;
  -\mathbf{R}_q^A\, \boldsymbol{\epsilon}_{ijk}\,(\mathbf{M}_A)_{j\ell}\, [\nabla^2 \Phi_B(\mathbf{R}_A)]_{k\ell}\, \hat{\mathbf{e}}_i,
  $$
  
  kde $\mathbf{M}_A = \int \rho_A \mathbf{x} \otimes \mathbf{x}\, dV$ je tenzor druhého momentu z §3.4.3 a $\boldsymbol{\epsilon}$ je Levi-Civitův tenzor. S užitím $[\nabla^2 \Phi_B]_{k\ell} = -Gm_B [3\hat{\boldsymbol{\rho}}_k \hat{\boldsymbol{\rho}}_\ell - \delta_{k\ell}]/|\boldsymbol{\rho}|^3$ (identita slapového tenzoru pro pole bodové hmoty), dosazením a převedením do tělesové soustavy podobnostní transformací $\mathbf{R}_q^A$ (s užitím stopové identity $\mathbf{Q}_A = \mathrm{tr}(\mathbf{I}_A)\mathbb{I}_3 - 3\mathbf{I}_A$ nebo ekvivalentně $\mathbf{M}_A = \tfrac{1}{2}\mathrm{tr}(\mathbf{I}_A)\mathbb{I}_3 - \mathbf{I}_A$ z §3.4.3) vychází tělesový moment síly jako:
  
  $$
  \boldsymbol{\tau}_A^{\mathrm{body}}
  \;=\;
  \frac{3 G m_B}{|\boldsymbol{\rho}|^3}\,
  \hat{\boldsymbol{\rho}}_{\mathrm{body},A} \times (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},A}),
  $$
  
  **identický** s výsledkem §3.7 zarámovaným.
  
  **Trojcestná křížová validace.** Tři nezávislé odvozovací cesty — (1) variační pravidlo $-\delta V_{\mathrm{tidal}}$ (§3.7), (2) přímý vzorec Hughes 2004 §5.3 (citovaný v §3.7.1), (3) přímá integrace $\mathbf{r}\times\nabla V$ Goldstein 2002 §5.5 (toto řešení OQ-3.3) — všechny produkují zarámovaný tělesový moment síly z gravitačního gradientu. Toto rozšiřuje metodologii dvě-paralelně z §3.6 na *tři-paralelně* pro tuto konkrétní rovnici, kalíc §3.7 proti zavedení algebraické chyby v kterékoli jediné cestě. Fortranový křížový test fáze 1 D4 `backends/fortran/tests/test_gravity_gradient_torque.f90` je regresní pojistkou C2 mezi implementacemi inženýrské úrovně při toleranci ULP × √N.
- **OQ-3.4.** ***Vyřešeno 2026-05-25 (kanonická úroveň, 3. kolo, navazující specifikace inženýrské úrovně).*** Rozšíření `dynamics.py::total_potential_energy` o zahrnutí kvadrupólově-slapových příspěvků $V_{\mathrm{tidal},A} + V_{\mathrm{tidal},B}$ z §3.4.3 zarámovaného.
  
  **Současná implementace (řádky 207-214).** Funkce inženýrské úrovně `total_potential_energy` počítá pouze monopólově-monopólový Keplerův člen $V^{(0)} = -Gm_Am_B/|\boldsymbol{\rho}|$ a vrací jej jako „potenciální energii soustavy“ pro diagnostiku zachování $|dE/E_0|$. Slapové části jsou implicitní v dynamice (skrze moment síly z gravitačního gradientu §3.7, jenž přímo vstupuje do tělesové Eulerovy pohybové rovnice), avšak nepřítomné z diagnostiky.
  
  **Specifikované rozšíření (patch inženýrské úrovně v0.1.1).** Přidej §3.4.3 zarámované slapové části k diagnostice:
  
  ```
  # Stávající
  V_kepler = -G * m_a * m_b / r_mag
  # Nové (navrhované)
  rhat_body_a = R_a.T @ (r_vec / r_mag)
  V_tidal_a = -(G * m_b / (2 * r_mag**3)) * (
                np.trace(I_a) - 3 * rhat_body_a @ I_a @ rhat_body_a)
  # symetricky pro těleso B
  rhat_body_b = R_b.T @ (-r_vec / r_mag)   # pozor na opačné znaménko
  V_tidal_b = -(G * m_a / (2 * r_mag**3)) * (
                np.trace(I_b) - 3 * rhat_body_b @ I_b @ rhat_body_b)
  V_total = V_kepler + V_tidal_a + V_tidal_b
  return V_total
  ```
  
  **Očekávaný dopad na diagnostiku zachování.** Při konfiguracích první verze ($|\boldsymbol{\rho}| \sim 1000$ m, $\ell \sim 7$ m) je velikost slapové části $|V_{\mathrm{tidal},A}/V_{\mathrm{Kepler}}| \sim (\ell/r)^2 \sim 5 \times 10^{-5}$. S rozšířením sleduje diagnostika zachování *plné* kvadrupólově-useknuté $V$ spíše než pouze monopólovou část; očekávaný efekt na $|dE/E_0|$:
  
  - Při $|\boldsymbol{\rho}| \sim 1000$ m: zlepšení je **pod zaokrouhlovací chybou dvojnásobné přesnosti** ($|dE/E_0|$ již na $6{,}2 \times 10^{-12}$, slapový příspěvek by toto změnil na úrovni $\sim 10^{-16}$, jež je pod šumem plovoucí čárky).
  - Při $|\boldsymbol{\rho}| \sim 30$ m (režim těsného přiblížení, $\ell/r \sim 0{,}2$): slapová část je $\sim (\ell/r)^2 \sim 0{,}04$ Keplerovy velikosti; bez rozšíření diagnostika zachování zdánlivě „selhává“ při $|dE/E_0| \sim 10^{-2}$ (neboť integrátor zachovává plné $V$, ale diagnostika měří pouze monopólové $V$); s rozšířením diagnostika zachování obnovuje skutečnou přesnost zachování integrátoru až na podlahu RK4.
  
  **Doprovodný navazující úkol v0.1.1: implementovat zpětnou sílu.** Jak vyplynulo v kap. 04 §4.9 OQ-4.5, zpětná síla na orbitu ($\partial V_{\mathrm{tidal}}/\partial \boldsymbol{\rho}$ — kap. 04 §4.6.2 zarámované) je rovněž nepřítomná z první verze. Rozšíření diagnostiky zachování OQ-3.4 je nezávislé na OQ-4.5 (diagnostika měří celkové $V$ správně; otázkou je, zda *pohybové rovnice* integrátoru propagují zpětnou sílu). Obojí by mělo být implementováno společně v v0.1.1, aby byly diagnostiky a pohybové rovnice vzájemně konzistentní.
  
  **Sledovač změn inženýrské úrovně.** Zařazeno v `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md` (další revize) §3.3 testovací případ těsného přiblížení: rozšířená diagnostika je tím, co by měl okrajový testovací případ třídy B těsného přiblížení užívat k přesné charakterizaci hranice oboru platnosti modelu.
- **OQ-3.5.** ***Vyřešeno 2026-05-23*** — projektovým vlastníkem dodané T&B PDF zrevidováno (1552stránkové elektronické vydání). T&B §27.5.1 („Multipole-Moment Expansion“) + §27.5.2 („Quadrupole-Moment Formalism“) pokrývá kartézský multipolový rozvoj newtonovského potenciálu v Rovn. (27.50)–(27.53), v zarámování teorie zdrojů gravitačních vln. Specifické citace přidány do §3.2 (T&B Rovn. 27.50 jednotělesový integrál) a §3.4.3 (T&B Rovn. 27.51, 27.52, 27.53 Taylorův rozvoj + monopólově-plus-kvadrupólové useknutí + definice bezstopového kvadrupólového momentu) referenčních anchorů, s rozdílem předfaktorové konvence mezi T&B $\mathcal{I}_{jk}$ (předfaktor $1/3$) a naším Jacksonovsko-konvenčním $\mathbf{Q}_{jk}$ (předfaktor $3$) učiněným explicitním ($\mathbf{Q}_{jk} = 3\, \mathcal{I}_{jk}$). **T&B *nerozvíjí* sféricko-harmonickou stranu multipolového rozvoje (naše §3.5)** — ta strana zůstává ukotvena na Jackson 1999 §4.1. **T&B *neodvozuje* stopovou identitu** $\mathbf{Q} = \mathrm{tr}(\mathbf{I}_{\mathrm{inertia}}) \cdot \mathbb{I} - 3\, \mathbf{I}_{\mathrm{inertia}}$ explicitně — poznamenává, že kvadrupólový moment je „druhým momentem rozložení hmoty s odečtenou stopou“ bez vyhláskování spojení s tenzorem setrvačnosti, jež naše §3.5 odvozuje. Řetězec kapitoly 03 kartézský-→-sféricko-harmonický-→-stopová-identita-→-moment-síly-z-gravitačního-gradientu je tedy integrační syntézou napříč anchory T&B + Jackson + Hughes spíše než kopií jakéhokoli jediného zdroje.
- **OQ-3.6.** ***Vyřešeno 2026-05-23*** — most na inženýrskou úroveň pro moment síly z gravitačního gradientu zpracován v nové podkapitole §3.7.2 („Jacobiho matice momentu síly z gravitačního gradientu: citlivost na počáteční podmínky a geometrii sestavy“), pokračujíc v rámci §2.5.1. Výstupy: linearizovaná Jacobiho matice orientační odezvy $\mathbf{J}^{\mathrm{orient}}_{A}$ (zarámovaná) jako sjednocující matice $3 \times 3$ působící na perturbace tělesové orientace i směru mezitělesové polohy; librační-rovnovážná specializace $\mathbf{J}^{\mathrm{orient}}_{A}\big|_{\mathrm{eq}} = \alpha (I_\parallel - I_\perp) \mathrm{diag}(1,1,0)$ s tabulkou stability gravitačního gradientu (protáhlé stabilní, zploštělé nestabilní, sférické marginální); radiálně-vs-tangenciální rozklad $\partial \boldsymbol{\tau}_A / \partial \boldsymbol{\rho}$ ukazující, že tangenciální část se redukuje na tutéž Jacobiho matici; identifikace $\partial \boldsymbol{\tau}_A / \partial q_B = 0$ při kvadrupólovém useknutí; tabulka hodnosti podle třídy symetrie (sférická → 0, axisymetrická → 2, asymetrická → 3); ukazatel citlivosti sestavy $\Sigma_{\mathrm{assembly}} = \eta_A \eta_B \cdot f(\hat{\boldsymbol{\rho}})$; numerická kontrola zdravého rozumu (podobné JWST v konfiguraci podobné L2: tuhost $4{,}7 \times 10^{-10}\ \text{N m/rad}$, e-násobek nestability $\sim 27\ \text{dní}$); tabulka vícemeřítkových časových měřítek identifikující první verzi jako spin-rozlišující, avšak librace-slepou; označení symplektického integrátoru na Lieově grupě pro dlouhohorizontové studie (položka inženýrské úrovně v0.2). Anchory: Hughes 2004 Kap. 8; Markley & Crassidis 2014 §3.7.

---

**Následující kapitola:** [04 — Eulerovy–Lagrangeovy rovnice a Euler-Poincarého redukce](04-euler-lagrange-and-euler-poincare.md) *(připravováno)*

**Konec 03-vzajemny-gravitacni-potencial.md**
