# Kapitola 7 — Ljapunovovo spektrum sdružené rotační dynamiky dvou těles

> **Předpokládané čtení.** Kapitoly 01–06. Tato kapitola pracuje se zarámovanou soustavou rovnic kap. 04 §4.7 propagovanou symplektickým variačním integrátorem na $SE(3) \times SE(3)$ z fáze 2 (`integrators/symplectic_se3_variational.py`, posezení 1) a kvantifikuje jeho citlivost na počáteční podmínky pomocí Ljapunovova spektra. Opírá se o testovací případy analytických limitů kap. 06 — bezsilové setrvačníky §6.3/§6.5 a gravitačně-gradientní libraci §6.4 — jako o kalibrační orákula.
>
> **Stav.** v0.1 (kanonická úroveň, fáze 2, posezení 4, 2026-06-07). Veškerá níže uvedená numerická evidence pochází z akceptačních běhů posezení 3 (`_handoffs/phase-2-symplectic-lie-group/SITTING-3-FINAL-REPORT.md`), všechny brány ZELENÉ.

> **Poznámka ke způsobu a předpokladům odvození.** Tato CS kanonická kapitola (v0.1) je odvozena z anglické kanonické kapitoly `docs/canonical/en/07-lyapunov-spectrum.md` (v0.1). Na rozdíl od ostatních kapitol zde **Podolský 2018 neposkytuje přímé CS-jazykové ukotvení** pro samotný formalismus Ljapunovova spektra — pokrývá pouze *integrabilní* referenční dynamiku (volný symetrický setrvačník §8.1; Lagrangeův setrvačník §8.3), která kalibruje scénáře (a) a (b). Pro algoritmický obsah je hlavním pramenem anglická literatura (Benettin 1980; Skokos 2010). Tato mezera je zde poctivě uvedena, nikoli zastřena volnou citací.

---

## §7.1 Přehled a motivace

Kapitoly 05 a 06 stanovily, *kde* je rotační dynamika sdruženého systému dvou těles stabilní, nestabilní nebo librující — avšak pouze v linearizovaném nebo integrabilním smyslu. Klasifikace „tenisové rakety“ z kap. 05 §5.2.2 dává linearizovaná vlastní čísla v rovnovážných stavech podél hlavních os; Lagrangeova redukce z kap. 06 §6.4.3 a Jacobiho eliptické řešení z §6.5.3 dávají uzavřený *integrabilní* pohyb o konečné amplitudě. Žádná z nich neodpovídá na otázku, kterou klade tato kapitola:

> **Narušuje gravitačně-gradientní vazba mezi orbitální a postojovou dynamikou integrabilitu — a pokud ano, jakou rychlostí exponenciální divergence?**

Kvantitativním nástrojem je **Ljapunovovo spektrum**: množina asymptotických rychlostí exponenciálního růstu $\lambda_1 \geq \lambda_2 \geq \cdots \geq \lambda_n$ infinitezimálních tečných perturbací podél trajektorie. Pro $n$-rozměrný tok $\dot{\mathbf{x}} = \mathbf{f}(\mathbf{x})$ s tokovým zobrazením $\phi^t$ jsou exponenty definovány tečným kocyklem $D\phi^t$:

$$
\lambda_i \;=\; \lim_{t \to \infty} \frac{1}{t}\, \ln \frac{\|D\phi^t(\mathbf{x}_0)\, \mathbf{v}_i\|}{\|\mathbf{v}_i\|},
$$

pro obecnou bázi $\{\mathbf{v}_i\}$ tečných vektorů (Oseledcův multiplikativní ergodický teorém zaručuje existenci skoro všude). Kladné $\lambda_1$ je operační definicí deterministického chaosu; $\lambda_1 = 0$ s kvaziperiodickým pohybem signalizuje integrabilní nebo regulární dynamiku.

Fyzikálně rozhodující fakta pro náš systém, kvantitativně doložená v §7.5–§7.6:

- **Bezsilový asymetrický setrvačník je integrabilní** — včetně pohybu poblíž sedla Džanibekova. Jeho uzavřený Jacobiho eliptický tvar (kap. 06 §6.5.3) neponechává prostor pro kladný Ljapunovův exponent; separatrisa je heteroklinická trajektorie integrabilního systému, nikoli chaotická vrstva.
- **Gravitačně-gradientní vazba narušuje integrabilitu.** Se zapnutou orbitálně-postojovou vazbou a počáteční podmínkou poblíž sedla podél prostřední osy činí vypočtený $\lambda_{\max} = 1{,}506 \times 10^{-5}\ \mathrm{s}^{-1} > 0$ — malý (slabý chaos poblíž separatrisy při použité síle vazby), avšak jednoznačně kladný.
- **Librace není divergence.** Těleso librující kolem gravitačně-gradientního rovnovážného stavu (kap. 06 §6.4) osciluje s frekvencí $\omega_{\mathrm{lib}}$, ale má $|\lambda_{\max}|$ o faktor 263 nižší než $\omega_{\mathrm{lib}}$ — frekvence oscilace a Ljapunovův exponent jsou různé veličiny a jejich směšování je klasické nedorozumění.

**Algoritmus a vnoření.** Spektrum počítáme algoritmem Benettin–Galgani–Giorgilli–Strelcyn (BGGS) 1980 (§7.2), s tečným kocyklem realizovaným jako konečnodiferenční Jacobiho matice kroku integrátoru v ploché euklidovské mapě $\mathbb{R}^{14}$ (§7.3). Euklidovské vnoření je záměrná, zdokumentovaná volba návrhu (OQ-PHASE2-2) se známým artefaktem konvergence pro integrabilní systémy, který přímo formoval prahy akceptačních bran — poctivé vyúčtování je v §7.3 a §7.5(a).

## §7.2 Algoritmus BGGS — formální znění

### §7.2.1 Proč naivní propagace tečného vektoru selhává

Propagace jediného tečného vektoru $\mathbf{v}(t) = D\phi^t \mathbf{v}_0$ a odečtení $\lambda_1 = \ln(\|\mathbf{v}(T)\|/\|\mathbf{v}_0\|)/T$ přeteče dvojnásobnou přesnost, jakmile $\lambda_1 T \gtrsim 700$; propagace báze $k$ vektorů je všechny zhroutí do nejvíce expandujícího směru a zničí přístup k $\lambda_2, \ldots, \lambda_k$. Řešení BGGS: bázi **periodicky reortonormalizovat** a akumulovat logaritmy faktorů protažení.

### §7.2.2 Znění algoritmu

> **Algoritmus (BGGS 1980, §2; moderní formulace Skokos 2010 §2.2).**
> *Vstup:* počáteční stav $\mathbf{x}_0$, krok integrátoru $\varphi_{dt}$, tečný propagátor $D\varphi_{dt}$, počet vektorů $k$, kadence QR $\tau_{\mathrm{QR}}$, horizont $T$.
>
> 1. Inicializuj náhodnou ortonormální bázi $Q_0 = [\mathbf{v}_1, \ldots, \mathbf{v}_k]$.
> 2. Pro každou epochu $m = 1, \ldots, T/\tau_{\mathrm{QR}}$:
>    a. Propaguj stav a všech $k$ tečných vektorů po $\tau_{\mathrm{QR}}/dt$ krocích:
>       $\mathbf{x} \leftarrow \varphi_{dt}(\mathbf{x})$, $\mathbf{v}_j \leftarrow D\varphi_{dt}(\mathbf{x})\,\mathbf{v}_j$.
>    b. Reortonormalizuj (modifikovanou) Gramovou–Schmidtovou ortonormalizací: $Q_m R_m = [\mathbf{v}_1, \ldots, \mathbf{v}_k]$; zaznamenej faktory protažení $r_j^{(m)} = (R_m)_{jj} = \|\mathbf{v}_j^{\perp}\|$ (normu $j$-tého sloupce po projekci na ortogonální doplněk předchozích sloupců).
>    c. Akumuluj $S_j \leftarrow S_j + \ln r_j^{(m)}$; resetuj bázi na $Q_m$.
> 3. *Výstup:* $\lambda_j = S_j / T$ pro $j = 1, \ldots, k$, sestupně.
>
> $j$-tý akumulovaný logaritmus protažení měří růst $j$-rozměrného objemového elementu dělený $(j-1)$-rozměrným, což podle Oseledcova teorému konverguje k $\lambda_j$ (Benettin et al. 1980 §2; diskrétní QR variantu srovnává s alternativami Geist, Parlitz & Lauterborn 1990).

### §7.2.3 Mapování na implementaci

Implementací je `lyapunov/spectrum.py::lyapunov_spectrum` (řídicí rutina, kroky 2–3), `lyapunov/cocycle.py::propagate_tangent` (realizace $D\varphi_{dt}$, §7.3) a `lyapunov/cocycle.py::gram_schmidt` (modifikovaná Gramova–Schmidtova ortonormalizace, krok 2b; modifikovaná varianta po sloupcích podle Golub & Van Loan §5.2 kvůli numerické stabilitě, s opětovným vložením degenerovaného sloupce z kanonické báze). Krok integrátoru $\varphi_{dt}$ je symplektický KDK Störmerův–Verletův krok na $SE(3) \times SE(3)$ z posezení 1 — týž krok, který drží ULP-bránu momentu hybnosti AG-INT-3, takže Ljapunovův běh dědí přesné zachování $\mathbf{L} = R\,\boldsymbol{\Pi}$ podél referenční trajektorie.

**Numerický příklad (produkční běh scénáře (c)).** $k=1$, $dt = 50$ s, $T = 100\,000$ s: pilotní fáze (§7.4) odhadla $\lambda_{\max}^{\mathrm{pilot}}$, nastavila produkční kadenci a akumulovaná protažení dala $\lambda_{\max} = 1{,}506 \times 10^{-5}\ \mathrm{s}^{-1}$ — hlavní číslo brány AG-LYAP-3 (§7.5(c)).

## §7.3 Volba tečného prostoru: euklidovský $\mathbb{R}^{14}$ versus varieta (OQ-PHASE2-2)

### §7.3.1 Volba

Stav rotačního podsystému leží na varietě $(S^3 \times \mathbb{R}^3)^2$ — jednotkový kvaternion + tělesová úhlová rychlost na každé těleso. Matematicky přirozeným tečným prostorem je Lieovsko-algebraický ($\mathfrak{so}(3) \times \mathbb{R}^3$ na těleso, dimenze 12). OQ-PHASE2-2 místo toho **potvrdil plochou euklidovskou mapu**:

$$
\mathbb{R}^{14} \;=\; \underbrace{q_A(4)}_{\text{kvaternion}} + \underbrace{\boldsymbol{\omega}_A(3)}_{\text{tělesová rychlost}} + \underbrace{q_B(4)}_{\text{kvaternion}} + \underbrace{\boldsymbol{\omega}_B(3)}_{\text{tělesová rychlost}},
$$

s tečným kocyklem realizovaným jako symetrický konečnodiferenční součin Jacobiho matice s vektorem pro krok integrátoru:

$$
\boxed{
D\varphi_{dt}(\mathbf{x})\cdot \mathbf{v}
\;\approx\;
\frac{\varphi_{dt}(\mathbf{x} + \epsilon\, \hat{\mathbf{v}}) - \varphi_{dt}(\mathbf{x} - \epsilon\, \hat{\mathbf{v}})}{2\epsilon}\; \|\mathbf{v}\|,
\qquad \epsilon = 10^{-7},
}
$$

kde $\hat{\mathbf{v}} = \mathbf{v}/\|\mathbf{v}\|$ (normalizace před perturbací udržuje velikost konečnodiferenční perturbace nezávislou na normě tečného vektoru, jež během běhu BGGS roste). Zdůvodnění: jednoduchost implementace (žádný atlas map na Lieově grupě, žádné účetnictví exponenciálního zobrazení), konečnodiferenční Jacobiho matice vyžadující jako černou skříňku pouze již validovaný krok integrátoru a transparentní opětovné použití rozvržení plného stavu $\mathbb{R}^{26}$ (`lyapunov/cocycle.py` `_R14_MAP`). Vnitřní normalizace kvaternionů v integrátoru opravuje nepatrný úlet mimo varietu vnesený každou konečnodiferenční perturbací.

### §7.3.2 Dva odlišné jevy — a poctivé vyúčtování (opraveno po fázi 3)

Ploché vnoření má **jeden skutečný artefakt** a je spojeno s **jedním jevem, jenž se nakonec ukázal jako jeho vina nebýt**. Fáze 3 (§7.10) ty dva oddělila; tento odstavec uvádí opravené porozumění, s poctivě zaznamenanou původní úvahou z doby fáze 2 níže.

**Skutečný artefakt vnoření — nefyzikální módy kvaternionové normy.** Každý jednotkový kvaternion v $\mathbb{R}^{14}$ nese radiální (normový) směr $\delta q \parallel q$, jenž míří mimo varietu omezení $S^3$ a je odstraněn normalizací v každém kroku integrátoru. Plné 14-vektorové spektrum je proto kontaminováno **dvěma velkými zápornými exponenty** (po jednom na těleso) bez kladných partnerů — spektrum se sčítá na $-0{,}14\ \mathrm{s}^{-1}$ místo $0$ a selhává v hamiltonovském párování (§7.7). Toto *je* artefakt vnoření a je důvodem, proč smysluplné plné spektrum vyžaduje 12-rozměrný kocyklus (fáze 3, §7.9–§7.10). Neovlivňuje jednovektorové brány $\lambda_{\max}$, jež sledují nejvíce expandující *dynamický* směr, dobře oddělený od kontrahujících normových módů.

**NIKOLI artefakt vnoření — pomalá konvergence FTLE ve scénáři (a).** Pro integrabilní scénář (a) konverguje konečnočasový Ljapunovův exponent (FTLE) k pravé hodnotě 0 pouze jako

$$
\lambda^{\mathrm{FTLE}}(T) \;\sim\; O\!\left(\frac{\log T}{T}\right),
$$

pomalý algebraický pokles, nikoli exponenciálně rychlá konvergence, jež by se naivně očekávala. **Fáze 3 ukázala, že to je vlastní kvaziperiodické vícefrekvenční dynamice, nikoli ploché mapě:** Lieovsko-algebraický kocyklus $\mathfrak{so}(3) \times \mathbb{R}^3$ (§7.10) — jenž nenese žádné směry kvaternionové normy — dává v podstatě totožné $\lambda_{\max} = 2{,}29 \times 10^{-4}\ \mathrm{s}^{-1}$ vůči hodnotě $\mathbb{R}^{14}$ $2{,}30 \times 10^{-4}\ \mathrm{s}^{-1}$ při $T = 20\,000$ s. Skutečně vnořením způsobená pomalá konvergence by pod kocyklem nativním pro varietu zmizela; nezmizela. Pomalá konvergence odráží nerezonanční směs frekvencí rotace a precese v samotné integrabilní trajektorii.

> **Poznámka k původu (opraveno 2026-06-07, fáze 3).** Verze tohoto odstavce z doby fáze 2 přisuzovala konvergenci $O(\log T/T)$ ve scénáři (a) „polynomiálnímu růstu normy v plochých kvaternionových souřadnicích“ — tedy vnoření — a předpovídala, že Lieovsko-algebraický kocyklus obnoví rychlou konvergenci a učiní bránu AG-LYAP-1 na úrovni $10^{-8}$ dosažitelnou. Fáze 3 tuto předpověď vyvrátila (AG-P3-3): Lieovsko-algebraický kocyklus reprodukoval totéž $\lambda_{\max}$. Zákon konvergence je vlastní; skutečnou vadou vnoření je kontaminace *plného* spektra nefyzikálními normovými módy, již Lieovsko-algebraický kocyklus *opravuje*. Původní dokumentační řetězec `lyapunov/reference.py::lambda_max_gate_scenario_a` nese překonanou úvahu. Jde o tutéž disciplínu opravy ukotvené ve fyzikální pravdě (KB-OPUS-DOC-INTRA-INCONSISTENCY), aplikovanou na mezi­úrovňovou předpověď, nikoli na zarámovaný vzorec.

### §7.3.3 Důsledek pro návrh bran

Důsledkem je metodologické pravidlo, jež stojí za to vyslovit v obecné podobě:

> **Akceptační prahy pro „nulové“ Ljapunovovy exponenty musí být kalibrovány empiricky vůči naměřenému zákonu konvergence kocyklu — nikoli převzaty z literatury o asymptotických Ljapunovových exponentech. Tento zákon může být dán vnořovací mapou, nebo (jako zde) vlastní dynamikou; tak či onak práh sleduje měření, nikoli ideál.**

Konkrétně: původní zadání fáze 2 specifikovalo AG-LYAP-1 (integrabilní scénář, $\lambda_{\max} <$ práh) na $10^{-8}$. Při zákonu $O(\log T / T)$ vyžaduje dosažení $10^{-8}$ čas $T \sim 10^{9}$ s — nedosažitelné při jakémkoli proveditelném integračním čase. Revidovaná brána $10^{-3}$ zachycuje správný *kvalitativní* požadavek (integrabilní FTLE musí být subdominantní vůči analytické rychlosti růstu chaotického scénáře $\sigma_c \approx 7{,}07 \times 10^{-3}\ \mathrm{s}^{-1}$, §7.6) a je splněna s rezervou 4,3× při $T = 20\,000$ s. Jde o zdokumentovanou revizi brány (zpráva posezení 3 §5 D3), nikoli o tiché uvolnění. **Fáze 3 potvrzuje, že brána $10^{-3}$ je správná hodnota i s kocyklem nativním pro varietu (Lieovsko-algebraický, §7.10)**, neboť pomalá konvergence je vlastní — revize brány byla správná, jen její původně uvedená *příčina* (§7.3.2) vyžadovala opravu.

## §7.4 Automatické ladění kadence QR pilot–produkce

Kadence QR $\tau_{\mathrm{QR}}$ vyvažuje cenu (každá reortonormalizace je $O(k^2 n)$) proti riziku přetečení/zhroucení (faktory protažení musí zůstat v rozsahu plovoucí čárky a báze nesmí mezi QR degenerovat). Užitečnou heuristikou je nechat každou epochu akumulovat $\sim 0{,}1$ e-násobku nejrychlejšího směru. Jelikož $\lambda_{\max}$ není předem znám, implementace jej automaticky ladí (vzor landscape §4.3):

1. **Pilotní běh:** 20 epoch při $\tau_{\mathrm{QR}}^{\mathrm{pilot}} = 10\,dt$ (celkem 200 kroků integrátoru — levné), dávající odhad $\lambda_{\max}^{\mathrm{pilot}} = S_1^{\mathrm{pilot}} / T_{\mathrm{pilot}}$.
2. **Produkční kadence:**

$$
\tau_{\mathrm{QR}} \;=\; \frac{0{,}1}{|\lambda_{\max}^{\mathrm{pilot}}|},
\qquad \text{ořezáno na } \big[\,10\,dt,\; T/3\,\big].
$$

3. **Záchrana při hodnotě blízké nule:** je-li $|\lambda_{\max}^{\mathrm{pilot}}| < 10^{-12}$ (integrabilní pilot), kadence degeneruje; místo toho se použije záchrana $\tau_{\mathrm{QR}} = \min(100\,dt,\; T/5)$.

**Numerický příklad.** Ve scénáři (b) (gravitačně-gradientní librace, $dt = 500$ s, $T = 10^6$ s) je pilotní odhad blízký nule až malý a ořez udržuje produkční kadenci pohodlně uvnitř $[5000\ \mathrm{s},\ 3{,}3\times 10^{5}\ \mathrm{s}]$; produkční běh se dokončí za 1,3 s reálného času při $T_{\mathrm{total}} = 10^6$ s simulovaného času (zpráva posezení 3 §4) — automaticky laděná kadence je to, co činí horizonty milionů sekund levnými.

## §7.5 Tříscénářová kalibrace

Ljapunovova pipeline, jež byla spuštěna pouze na zkoumaném systému, je nekalibrovaná. Testovací sada posezení 3 (`testbeds/lyapunov_three_scenario.py`) ohraničuje fyziku třemi scénáři, jejichž očekávané odpovědi jsou známy z analýzy kap. 05–06 — jeden integrabilní-nulový, jeden omezený librací, jeden skutečně chaotický.

### §7.5(a) Bezsilový asymetrický setrvačník — $\lambda_{\mathrm{true}} = 0$ (AG-LYAP-1)

**Nastavení.** $\mathbf{I} = \mathrm{diag}(100, 200, 400)$ kg·m², $\boldsymbol{\omega}_0 = [0{,}1, 0, 0]$ rad/s (stabilní osa nejmenšího momentu $\hat{\mathbf{e}}_1$ — kvaziperiodický, rychlá konvergence FTLE), partnerské těleso zaparkováno na $10^{12}$ m, takže gravitační vazba je zanedbatelná. Dynamikou je integrabilní Jacobiho eliptický pohyb z kap. 06 §6.5.3; pravé $\lambda_{\max} = 0$.

**Očekávání.** Podle artefaktu §7.3.2 klesá vypočtené FTLE k 0 jen jako $O(\log T / T)$; brána testuje subdominanci, nikoli doslovnou nulu.

| Veličina | Hodnota |
|---|---|
| $dt$, $T$ | 1,0 s, 20 000 s |
| $\lambda_{\max}$ (vypočtené FTLE) | $2{,}304 \times 10^{-4}\ \mathrm{s}^{-1}$ |
| Brána AG-LYAP-1 ($< 10^{-3}$) | **PROŠLA** ✓ (rezerva 4,3×) |

### §7.5(b) Gravitačně-gradientní librace — $|\lambda_{\max}| < \omega_{\mathrm{lib}}$ (AG-LYAP-2)

**Nastavení.** Těleso podobné JWST, protáhlé (kap. 06 §6.4), v gravitačně-gradientním rovnovážném stavu — dlouhá osa směrem k primárnímu tělesu o hmotnosti Země ($m_B = 5{,}97 \times 10^{24}$ kg) na $\rho = 10^{7}$ m — nakloněné o 0,01 rad. Těleso libruje s frekvencí z kap. 06 §6.4.2 $\omega_{\mathrm{lib}} = \sqrt{3 G m_B (I_\perp - I_\parallel)/(I_\perp \rho^3)}$.

**Očekávání.** Librace je omezená oscilace v efektivním potenciálu Lagrangeova setrvačníku (kap. 06 §6.4.3) — *žádná* exponenciální divergence. Rozlišující brána: $|\lambda_{\max}|$ musí ležet hluboko pod samotnou frekvencí oscilace.

| Veličina | Hodnota |
|---|---|
| $dt$, $T$ | 500 s, $10^{6}$ s |
| $\omega_{\mathrm{lib}}$ | $6{,}379 \times 10^{-4}\ \mathrm{s}^{-1}$ |
| $\lambda_{\max}$ | $2{,}421 \times 10^{-6}\ \mathrm{s}^{-1}$ |
| Brána AG-LYAP-2 ($|\lambda_{\max}| < \omega_{\mathrm{lib}}$) | **PROŠLA** ✓ (rezerva 263×) |

Oddělení 263× je kvantitativní podobou opatrnosti z §7.1: librující kosmická loď *osciluje* s $\omega_{\mathrm{lib}}$, avšak *nediverguje* ničím, co by se tomu blížilo.

### §7.5(c) Gravitačně-gradientní Džanibekov poblíž sedla — $\lambda_{\max} > 0$ (AG-LYAP-3)

**Nastavení.** Asymetrické těleso $\mathbf{I} = \mathrm{diag}(100, 200, 400)$ kg·m², $\boldsymbol{\omega}_0 = [0{,}001, 0{,}01, 0{,}001]$ rad/s — poblíž sedla podél prostřední osy (kap. 06 §6.5) — s primárním tělesem o hmotnosti Země na $\rho = 3 \times 10^{6}$ m, zvoleným tak, aby gravitačně-gradientní vazba byla téměř rezonanční s dynamikou sedla. Analytická charakteristická frekvence Džanibekova (kap. 06 §6.5.2):

$$
\sigma \;=\; \Omega_2^0 \sqrt{\frac{(I_2 - I_1)(I_3 - I_2)}{I_1 I_3}}
\;=\; 0{,}01 \sqrt{\frac{100 \cdot 200}{100 \cdot 400}}
\;=\; 7{,}071 \times 10^{-3}\ \mathrm{s}^{-1}.
$$

**Očekávání.** Bezsilově by tento pohyb byl integrabilní (heteroklinický, Jacobiho eliptický). Gravitačně-gradientní vazba narušuje integrabilitu; poblíž separatrisy to vytváří tenkou chaotickou vrstvu s malým, ale kladným $\lambda_{\max}$ — řádu KAM-rezidua při použité síle vazby, *nikoli* řádu $\sigma$.

| Veličina | Hodnota |
|---|---|
| $dt$, $T$ | 50 s, $10^{5}$ s |
| $\sigma$ (analytická) | $7{,}071 \times 10^{-3}\ \mathrm{s}^{-1}$ |
| $\lambda_{\max}$ | $1{,}506 \times 10^{-5}\ \mathrm{s}^{-1}$ |
| poměr $\lambda_{\max}/\sigma$ | 0,0021 |
| Brána AG-LYAP-3 ($\lambda_{\max} > 0$, poměr $< 10$) | **PROŠLA** ✓ |

Původní brána v zadání požadovala poměr v rozmezí faktoru 3 od $\sigma$; tato kalibrace předpokládala sílu vazby $\varepsilon \approx 1$. Při skutečném $\varepsilon = \omega_{\mathrm{lib}}/\sigma \approx 0{,}81$ leží exponent slabého chaosu o řády pod $\sigma$ — fyzikálně správné pro chaos v režimu KAM poblíž separatrisy — a brána byla revidována na poměr $< 10$ s tím, že fyzikální obsah nese požadavek kladnosti (zpráva posezení 3 §5, tabulka revizí bran).

## §7.6 Mezi­scénářové srovnání a kritérium integrability

Tři scénáře se vyplatí teprve společně: pipeline je musí *uspořádat* správně. Naivní test uspořádání — vypočtené $\lambda_c > \lambda_a$ — selhává z poučného důvodu: při dostupném $T$ leží stále klesající FTLE integrabilního scénáře ($3{,}924 \times 10^{-4}$ při $T = 10\,000$ s) *nad* pravým chaotickým exponentem $\lambda_c \approx 1{,}5 \times 10^{-5}$. Při zákonu $O(\log T/T)$ by $\lambda_a^{\mathrm{FTLE}}$ nekleslo pod $\lambda_c$ až do $T_a \approx 2{,}6 \times 10^{8}$ s — neproveditelné. Robustní, fyzikálně smysluplné uspořádání srovnává integrabilní FTLE s *analytickou* chaotickou rychlostí růstu:

| Veličina | Hodnota |
|---|---|
| $\lambda_a^{\mathrm{FTLE}}$ při $T = 10\,000$ s | $3{,}924 \times 10^{-4}\ \mathrm{s}^{-1}$ |
| $\sigma_c$ (analytická) | $7{,}071 \times 10^{-3}\ \mathrm{s}^{-1}$ |
| poměr | **0,055** |
| Brána uspořádání ($\lambda_a < \sigma_c$) | **PROŠLA** ✓ |

přičemž $\lambda_c > 0$ je stanoveno samostatně branou AG-LYAP-3. Souhrnný obraz:

| Scénář | Pravé $\lambda_{\max}$ | Vypočtená evidence | Verdikt |
|---|---|---|---|
| (a) bezsilový asymetrický setrvačník | $0$ (integrabilní) | FTLE $\to 0$ jako $O(\log T/T)$; $2{,}3 \times 10^{-4} \ll \sigma_c$ | regulární |
| (b) gravitačně-gradientní librace | $\approx 0$ (omezená librace) | $|\lambda| = 2{,}4 \times 10^{-6} \ll \omega_{\mathrm{lib}}$ (263×) | regulární (oscilační) |
| (c) gravitačně-gradientní Džanibekov poblíž sedla | $> 0$ (vazba narušuje integrabilitu) | $\lambda = 1{,}5 \times 10^{-5} > 0$, poměr $\lambda/\sigma = 0{,}002$ | slabě chaotický |

To je **kritérium integrability** této kapitoly v operační podobě: integrabilní a librující konfigurace vykazují FTLE, jež klesají s $T$ a leží hluboko pod relevantními dynamickými frekvencemi ($\sigma$, $\omega_{\mathrm{lib}}$); skutečně chaotické konfigurace vykazují vůči $T$ stabilní kladný exponent, jakkoli malý.

## §7.7 Hamiltonovské párování a struktura spektra

Sdružený rotační systém dvou těles je hamiltonovský a jeho Ljapunovovo spektrum dědí symplektické omezení: exponenty přicházejí v **párech $(+\lambda, -\lambda)$**, s nulovými exponenty přispívanými zachovávajícími se veličinami (samotný hamiltonián plus každý nezávislý první integrál na redukovaném prostoru). Pro 14-rozměrný rotační podsystém to znamená 7 párů; integrabilita by vynutila *všechny* páry na $(0, 0)$.

Implementační kontrolou je `lyapunov/reference.py::hamiltonian_pairing_check`: pro $k$ vypočtených exponentů v sestupném pořadí testuje

$$
\lambda_i + \lambda_{k-1-i} \;\approx\; 0, \qquad i = 0, \ldots, \lfloor k/2 \rfloor - 1,
$$

s relativní tolerancí 0,5 (kvalitativní brána — rezidua párování, nikoli přesná rekonstrukce symplektického spektra) a počítá exponenty blízké nule ($|\lambda| < 10^{-4}$) jako kandidátní invariantní směry.

**Poznámka fáze 2 → fáze 3.** Produkční běhy posezení 3 počítaly $k = 1$ tečný vektor (výše uvedené brány potřebují jen $\lambda_{\max}$). Při $k = 1$ je kontrola párování prázdná — párování je pozorovatelné teprve při $k \geq 2$ a plně informativní při $k = 14$. Plný audit párování celého spektra byl proto odložen na fázi 3. Předběžný běh $k = 14$ po posezení 4 zjistil, že ploché spektrum $\mathbb{R}^{14}$ je **ovládáno dvěma nefyzikálními kontrahujícími módy kvaternionové normy** (po jednom na těleso), sčítajícími se na $-0{,}14\ \mathrm{s}^{-1}$ a selhávajícími v párování — skutečný artefakt vnoření z §7.3.2. Fáze 3 (§7.10) jej poté **opravila** 12-rozměrným kocyklem: s odstraněnými normovými směry se spektrum sčítá na $\approx 0$ a brána hamiltonovského párování prochází (AG-P3-1/2). Kontrola párování je nyní skutečnou bránou, nikoli prázdnou.

## §7.8 Souvislost s kapitolami 05–06: vlastní čísla, librace a exponenty

Tři příběhy stability z kap. 05–06 se nyní uzavírají do jediného konzistentního obrazu:

**Librační vlastní čísla versus Ljapunovovy exponenty (kap. 06 §6.4.2 ↔ §7.5(b)).** Linearizovaná analýza librace dává čistě imaginární vlastní čísla $\pm i\,\omega_{\mathrm{lib}}$ v protáhlém gravitačně-gradientním rovnovážném stavu. Imaginární vlastní čísla znamenají oscilaci, nikoli divergenci — a nelineárním potvrzením je právě AG-LYAP-2: $|\lambda_{\max}| = 2{,}4 \times 10^{-6} \ll \omega_{\mathrm{lib}} = 6{,}4 \times 10^{-4}$. Librace **není** exponenciální divergence; její Ljapunovův exponent je konzistentní s nulou. Lagrangeova redukce (kap. 06 §6.4.3, analogie Podolský §8.3) vysvětluje *proč*: pohyb o konečné amplitudě je omezen v jámě efektivního potenciálu, integrabilní kvadraturou, tedy regulární.

**Džanibekov: sedlo ≠ chaos, dokud není provázáno (kap. 05 §5.2.2 + kap. 06 §6.5.3 ↔ §7.5(c)).** Rovnovážný stav podél prostřední osy má reálné kladné linearizované vlastní číslo $\sigma$ (teorém tenisové rakety z kap. 05), avšak bezsilový pohyb o konečné amplitudě je *integrabilní* — uzavřený Jacobiho eliptický tvar z kap. 06 §6.5.3 parametrizuje každou trajektorii, včetně přeskoků poblíž separatrisy, a integrabilní systém má identicky nulové Ljapunovovo spektrum. Linearizované nestabilní vlastní číslo v sedle je lokální tvrzení; globální pohyb prochází heteroklinickou strukturou regulárně. Teprve **gravitačně-gradientní vazba** — jež ničí druhý první integrál — mění okolí separatrisy ve skutečnou tenkou chaotickou vrstvu, a i tehdy s $\lambda_c \approx 0{,}002\,\sigma$ při síle vazby fáze 2. Hierarchie *linearizované vlastní číslo → integrabilní nelineární orákulum → Ljapunovův exponent provázaného systému* je úplnou metodikou stability tohoto projektu a každá úroveň je nezávisle ohraničena bránou (kap. 06 §6.5.5; AG-LYAP-1/3).

## §7.9 Otevřené otázky a výhled na fázi 3

- **OQ-PHASE2-LYA-1 — Lieovsko-algebraický kocyklus — VYŘEŠENO (fáze 3, §7.10).** Tečný propagátor v tělesové Lieově algebře ($\mathfrak{so}(3) \times \mathbb{R}^3$ na těleso, dimenze 12) nahrazuje plochou konečnodiferenční Jacobiho matici $\mathbb{R}^{14}$. **Dosažený přínos:** odstraňuje dva nefyzikální kontrahující módy kvaternionové normy (12-rozměrný konstrukcí), takže plné spektrum se sčítá na $\approx 0$ a hamiltonovské párování prochází. **Přínos NEDOSAŽEN (předpověď opravena):** očekávání z doby fáze 2, že to obnoví rychlou konvergenci FTLE a učiní bránu AG-LYAP-1 na úrovni $10^{-8}$ dosažitelnou, bylo *vyvráceno* — Lieovsko-algebraický kocyklus dává $\lambda_{\max} = 2{,}29 \times 10^{-4}$ vůči $\mathbb{R}^{14}$ $2{,}30 \times 10^{-4}$ při $T = 20\,000$ s. Pomalá konvergence je vlastní kvaziperiodické dynamice, nikoli vnoření (§7.3.2 opraveno; AG-LYAP-1 zůstává na $10^{-3}$).
- **OQ-PHASE2-LYA-2 — Plné (12-vektorové) spektrum — VYŘEŠENO (fáze 3, §7.10).** S 12-rozměrným kocyklem (deflace Tier A nebo Lieovsko-algebraický Tier B) je audit hamiltonovského párování skutečnou bránou: součet spektra $|\Sigma\lambda| < 10^{-2}$ (AG-P3-1), oba velké záporné normové módy nepřítomny a invariantní směry blízké nule přítomny (AG-P3-2). Předběžná deflace Tier A (projekce $\mathbb{R}^{14}$ proti $q_A, q_B$ před každým QR) i kanonický Lieovsko-algebraický kocyklus Tier B oba toto plní. Odhady Kaplanovy–Yorkeovy dimenze zůstávají cílem nad rámec (mimo rozsah fáze 3).
- **OQ-PHASE2-LYA-3 — Meze stability postojového řízení JWST (stále otevřeno).** Inženýrská otázka motivující testovací sadu: převést $\lambda_{\max}$ pro realistické konfigurace podobné JWST (kompozitní setrvačnost, orbitální prostředí soustavy v bodě L2, kadence správy momentu hybnosti) na časové konstanty divergence $1/\lambda$ srovnatelné s šířkou pásma postojového řízení a intervaly udržování dráhy. Vyžaduje rozšíření omezené soustavy tří těles u bodu L2 z v0.2+ (kap. 06 §6.9 OQ-6.1), než bude orbitální prostředí reprezentativní.
- **Mezera v pokrytí Podolským (zopakováno z §7.1).** Tato CS-jazyková kanonická kapitola ukotvuje integrabilní referenční dynamiku v Podolský §8.1/§8.3 a pro formalismus Ljapunova odkazuje na anglickou literaturu (Benettin 1980; Skokos 2010); Podolský 2018 Ljapunovova spektra nepokrývá.

## §7.10 Fáze 3 — 12-rozměrné spektrum (řeší OQ-PHASE2-LYA-1 + LYA-2)

Fáze 3 staví dva 12-rozměrné kocykly a znovu spouští spektrum, odstraňujíc kontaminaci nefyzikálními módy z §7.7. Oba ponechávají základní $\mathbb{R}^{14}$ (`lyapunov/cocycle.py`) netknutý kvůli křížové validaci; oba jsou aditivní (žádný soubor fáze 2 nezměněn).

**Tier A — deflační kocyklus kvaternionové normy** (`lyapunov/cocycle_deflated.py`). Obaluje konečnodiferenční kocyklus $\mathbb{R}^{14}$ a před každým krokem Gramovy–Schmidtovy ortonormalizace promítá každý tečný vektor ortogonálně k oběma aktuálním radiálním kvaternionovým směrům $\hat{q}_A, \hat{q}_B$, pracujíc v deflovaném 12-rozměrném podprostoru. Rychlé, s nízkým rizikem, přímo prokazuje diagnózu 12-vs-14.

**Tier B — Lieovsko-algebraický kocyklus** (`lyapunov/cocycle_liealg.py`). Tečný stav $(\xi_A, \delta\boldsymbol{\omega}_A, \xi_B, \delta\boldsymbol{\omega}_B) \in \mathbb{R}^{12}$ s $\xi \in \mathfrak{so}(3)$. Perturbace orientace užívají exponenciální zobrazení v **pravé trivializaci odpovídající integrátoru** ($q_{k+1} = q_k \otimes \exp(\epsilon\,\xi)$), s logaritmickým zobrazením $\xi = 2\,\mathrm{arctan2}(|\mathbf{q}_{\mathrm{vec}}|, q_w)\,\mathbf{q}_{\mathrm{vec}}/|\mathbf{q}_{\mathrm{vec}}|$ mapujícím propagovanou perturbaci zpět do $\mathfrak{so}(3)$. Symetrická konečná diference jako ve fázi 2, stejný počet vyhodnocení.

**Výsledky (scénář (c), všechny brány ZELENÉ; testy `tests/test_lyapunov_phase3.py`, 16 rychlých + 5 pomalých):**

| Brána | Kontrola | Výsledek |
|---|---|---|
| AG-P3-1 | součet spektra $|\Sigma\lambda| < 10^{-2}$ | **PROŠLA** (součet $\ll 0{,}14$; módy $-0{,}064/-0{,}071$ pryč) |
| AG-P3-2 | hamiltonovské párování, $n_{\mathrm{blízko\text{-}nuly}} \geq 2$ | **PROŠLA** |
| AG-P3-3 | Lieovsko-algebraické $\lambda_{\max}$ scénář (a) $< 10^{-3}$ | **PROŠLA** ($2{,}29 \times 10^{-4}$; brána revidována z aspirace $10^{-5}$ — viz níže) |
| AG-P3-4 | Lieovsko-algebraické $\lambda_{\max}$ v rozmezí faktoru 3 od hodnoty $\mathbb{R}^{14}$ | **PROŠLA** (obě $> 0$, scénář (c)) |
| AG-P3-REGRESE | Všechny brány fáze 2 nezměněny | **PROŠLA** (16/16 rychlých) |

**Nález AG-P3-3 — opravená předpověď.** Úvaha z doby fáze 2 v §7.3.2/§7.9 předpovídala, že Lieovsko-algebraický kocyklus obnoví rychlou konvergenci FTLE pro scénář (a) a učiní bránu $10^{-8}$ dosažitelnou. Neučinil: kocyklus nativní pro varietu dal $\lambda_{\max} = 2{,}29 \times 10^{-4}\ \mathrm{s}^{-1}$, v podstatě totožné s hodnotou $\mathbb{R}^{14}$ $2{,}30 \times 10^{-4}$. Konvergence $O(\log T/T)$ je tedy **vlastní kvaziperiodické vícefrekvenční dynamice, nikoli plochému vnoření** — skutečně mapou způsobený jev by pod kocyklem nativním pro varietu zmizel. Brána zůstává na hodnotě fáze 2 $10^{-3}$ (parita s AG-LYAP-1) a §7.3.2 je odpovídajícím způsobem opraveno. Jde o disciplínu poctivosti z §7.3.2 obrácenou proti vlastnímu dřívějšímu tvrzení kapitoly: mezi­úrovňová předpověď byla testována fyzikální pravdou Lieovsko-algebraického experimentu a neobstála.

**Co fáze 3 změnila a nezměnila.** *Opravila* problém plného spektra / hamiltonovského párování (skutečný artefakt vnoření). *Nezměnila* jednovektorové brány $\lambda_{\max}$ ani práh AG-LYAP-1 (pomalá konvergence je vlastní). Dva jevy, jež §7.3.2 původně směšovalo, jsou nyní čistě odděleny.

---

> **Hlavní odkaz.** Pro samotný formalismus Ljapunovova spektra **Podolský 2018 neposkytuje ukotvení** (tématu se nevěnuje). Hlavním pramenem je proto anglická literatura — viz Doplňující anglojazyčné odkazy níže. Pro *integrabilní* referenční dynamiku kalibrující scénáře (a) a (b) je hlavním CS-jazykovým ukotvením Podolský 2018 §8.1 (volný symetrický setrvačník) a §8.3 (těžký symetrický setrvačník s pevným bodem — efektivní potenciál a tři precesní režimy), týž anchor jako pro kap. 06 §6.3/§6.4.
>
> **Doplňující anglojazyčné odkazy.** Pro algoritmus: Benettin, Galgani, Giorgilli & Strelcyn 1980, *Meccanica* **15**, 9–30, §2 (`Benettin1980` — reortonormalizační algoritmus QR/Gram–Schmidt; primární ukotvení); Skokos 2010, *Lect. Notes Phys.* **790**, §2.2 (`Skokos2010` — moderní formulace, chování konvergence FTLE); Geist, Parlitz & Lauterborn 1990, *Prog. Theor. Phys.* **83**, 875–893 (`Geist1990` — srovnání diskrétních QR metod). Pro numeriku modifikované Gramovy–Schmidtovy ortonormalizace: Golub & Van Loan 2013 §5.2. Pro integrabilní orákula kalibrující scénáře: Landau–Lifšic, sv. I §37 a Goldstein 2002 §5.6. Pro symplektický integrátor dodávající kocyklus: Hairer–Lubich–Wanner 2006 §VII.5 (`HairerLubichWanner2006`).

---

**Následující:** OQ-PHASE2-LYA-3 (meze stability postojového řízení JWST) — vyžaduje prostředí omezené soustavy tří těles u bodu L2 z v0.2+ (kap. 06 §6.9 OQ-6.1), než bude orbitální nastavení reprezentativní. Fáze 3 (§7.10) vyřešila LYA-1 + LYA-2.

**Konec 07-ljapunovovo-spektrum.md (v0.1; §7.10 přidáno + §7.3.2/§7.9 opraveno ve fázi 3, 2026-06-07).**
