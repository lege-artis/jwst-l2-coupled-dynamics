# lege-artis/jwst-l2-coupled-dynamics — kanonická úroveň — Úvod (CS)

> **Publikum.** Matematik, fyzik, akademický recenzent nebo spolupracovník, který chce pochopit, *proč* kód pro sdruženou dynamiku v tomto projektu počítá to, co tvrdí, že počítá — odvozeno z prvních principů v moderní formě nezávislé na volbě souřadnicové soustavy.
>
> **Čas čtení.** Tento úvod přibližně 30 minut. Celá kanonická úroveň (tento + kapitoly 01–06) po dokončení přibližně 3–4 hodiny.
>
> **Licence.** CC-BY-SA-4.0 (záměr; projekt zatím není veřejně zrcadlen).
>
> **Doprovodné úrovně.** Inženýrská úroveň `docs/engineer/en/00-quick-start.md` (připravovaná) pokryje stejný rozsah z pohledu praktikanta používajícího knihovnu. Shad-úroveň `SHAD-NARRATIVE.md` (již odeslaná) pokrývá průvodce prvním prototypem psaný inženýrskou narací.

> **Poznámka ke způsobu a předpokladům odvození.** Tato CS kanonická úroveň (v0.1) je odvozena z anglické kanonické úrovně `docs/canonical/en/00-introduction.md` (v0.1). Terminologie je ukotvena v Podolský 2018 *Teoretická mechanika ve třech knihách* (Matfyzpress, ISBN 978-80-7378-499-7). Tělesový text je psán jako samostatný CS dokument; EN verze zůstává autoritativním pramenem při budoucích aktualizacích.

---

## §1. Rozsah kanonické úrovně

Tento projekt — dočasně uložený v adresáři `experiments/jwst-l2-first-cut/` v pracovním prostoru a zamýšlený jako základ budoucího repozitáře `lege-artis/jwst-l2-coupled-dynamics` — poskytuje referenční implementaci dynamiky dvou gravitačně provázaných tuhých těles. Primárním tělesem je těleso podobné JWST, sekundárním těleso sondy; systém se nachází v okolí (v konečném důsledku) bodu rovnováhy L2 soustavy Země–Slunce.

Kanonická úroveň odvozuje z prvních principů pohybové rovnice, které tento systém řídí. Odvozování probíhá v moderním jazyku nezávislém na volbě souřadnicové soustavy a ke konkrétním soustavám přechází teprve tehdy, kdy to výpočet vyžaduje.

**Rozsah v0.1 této kanonické úrovně** (rozsah zde zpracovaný):

| Prvek domény | Zpracování |
|---|---|
| Dvě tuhá tělesa (rozšířená, s plnými tenzory setrvačnosti) | Ano |
| Vzájemný gravitační potenciál mezi nimi | Ano, odvozeno z prvních principů (multipolový rozvoj) |
| Tuhá rotace kolem těžiště | Ano, odvozeno pomocí Lagrangeova formalismu + Euler-Poincarého redukce |
| Relativní orbitální pohyb v těžišťové soustavě | Ano |
| Dynamika u bodu L2 v omezené soustavě tří těles (Slunce–Země–kosmická loď) | **Odloženo na v0.2** této kanonické úrovně |
| Relativistické korekce 1PN/2PN | **Záměrně vyloučeny** (zdůvodnění v §6 níže; příslušný rychlostní poměr je $v/c \sim 10^{-5}$ pro orbitální měřítka JWST, příspěvek k rotační vazbě zanedbatelný) |
| Vnější perturbace (tlak slunečního záření, perturbace třetím tělesem Měsíce atd.) | **Mimo rozsah** kanonické úrovně; inženýrská úroveň je zkatalogizuje |
| Vnitřní pružnost / vibrační módy | **Mimo rozsah**; přiblížení tuhým tělesem je výchozí předpoklad kanonické úrovně |

**Soupis kapitol v0.1** (kapitoly, které tuto úroveň dokončují na v0.1):

| Kapitola | Název | Stav | CS-jazykové ukotvení (Podolský 2018) |
|---|---|---|---|
| 00 | Úvod (tento soubor) | **Zpracována** | — (úvodní) |
| 01 | Konfigurační varieta a variační princip | **Skeletována** | §6.2 + §6.4 (ortogonální transformační matice + mapa Eulerových úhlů jako elementární předchůdce; abstrakce Lieových grup je moderní Marsden-Ratiu, nikoli v Podolském) |
| 02 | Kinetická energie na $SE(3) \times SE(3)$ v intrinzické formě | **Zpracována** (v0.1, doplnění v 2. kole 2026-05-24 dop.: §2.5.1 citlivost d**I**/dt, §2.6.1 symplektická transformace těžiště, §2.7.1 levá invariance) | §6.3 Rovn. (6.4) $\boldsymbol{\Omega}$ tělesové báze + (6.8) zarámovaná identita pro časovou derivaci vektoru + (6.2)/(6.10) trivializace; §7.1 (7.3)/(7.7)/(7.9)/(7.11) tenzor setrvačnosti + (7.12) Steinerova věta |
| 03 | Vzájemný gravitační potenciál a multipolový rozvoj | **Zpracována** (v0.1, s opravou znaménka §3.7.2 OQ-FORT-1, 2026-05-24 odp.) | Částečně: §7.1 algebra tenzoru setrvačnosti jako most ke gravitačnímu kvadrupolovému tenzoru (§3.5); §8.3 librace těžkého symetrického setrvačníku s pevným bodem jako klasický precedent pro gravitačně-gradientní libraci (§3.7). Moderní materiál multipolového rozvoje a gravitačního gradientu je mimo rozsah Podolského. |
| 04 | Eulerovy–Lagrangeovy rovnice a Euler-Poincarého redukce | **Zpracována** (v0.1, 3. kolo kanonické úrovně 2026-05-25: §4.2 sestavení Lagrangiánu + §4.3 Přístup A Eulerovy–Lagrangeovy rovnice v mapě Eulerových úhlů + §4.4 Přístup B Euler-Poincarého redukce s explicitním odvozením $\mathrm{ad}^*$ + §4.5 křížová validace A≡B + §4.6 translační rovnice + zpětná síla + §4.7 zarámovaná soustava rovnic + řádková křížová validace vůči `dynamics.py`; všechny čtyři vnitřní OQ 4.1–4.4 vyřešeny; dopředu orientované OQ-4.5..4.8 pro v0.1.1/v0.2) | §7.2 zarámované (7.15) Eulerovy dynamické rovnice + §7.3 odvození plného Lagrangeova formalismu přes Eulerovy úhly + §6.4 (6.14)–(6.17) aparát mapy — přímé CS ukotvení pro Přístup A; metodologie dvou paralelních přístupů formalizovaná v této kapitole je implicitní v Podolského prezentaci §7.2 + §7.3 se dvěma odvozeními. Euler-Poincarého redukce (Přístup B) je moderní Marsden-Ratiu, nikoli v Podolském. |
| 05 | Eulerovy rovnice v tělesové souřadnicové soustavě jako specializace | **Zpracována** (v0.1, 3. kolo kanonické úrovně 2026-05-25: §5.2 specializace na hlavní osy s rozvojem po složkách + §5.3 odstavce pro třídy symetrie + §5.4 moment síly z gravitačního gradientu dosazený do tvaru s hlavními osami + §5.4.1 křížová kontrola libračního rovnovážného stavu vůči Ch 03 §3.7.2 + §5.5 zákony zachování pro případ bez momentu sil i s momentem sil + §5.6 ukazatele na testovací programy; všechny čtyři vnitřní OQ 5.1–5.4 vyřešeny; dopředu orientované OQ-5.5..5.6 pro v0.2) | §7.2 zarámované (7.15) — přímá CS analogie centrálního výsledku kapitoly; §8.1 (8.1)–(8.2) volný (bezsilový) symetrický setrvačník jako analytický testovací případ; §8.2 setrvačník se třením jako varianta s útlumem postoje |
| 06 | Analytické limity a validační brány | **Zpracována** (v0.1, 3. kolo kanonické úrovně 2026-05-25: §6.1 disciplína tříd A/B + §6.2 testovací případ limitu Keplera + §6.3 testovací případ volného symetrického setrvačníku (implementován v první verzi) + §6.4 testovací případ gravitačně-gradientní librace (třída A linearizovaná + třída A plně nelineární přes uzavřený vzorec) + §6.5 testovací případ asymetrického setrvačníku Džanibekova + §6.6 testovací případ oddělení při nekonečné vzdálenosti + §6.7 brána C2 bitové identity mezi implementacemi; dopředu orientované OQ-6.1..6.6 pro v0.2) | §8.1 volný setrvačník + §8.2 setrvačník se třením + §8.3 Lagrangeův setrvačník — všechny tři jsou analytické testovací případy s Podolského uzavřenými řešeními, přímo použitelnými jako orákula pro validační brány. |

Rozšíření **v0.2** přidá vrstvu omezené soustavy tří těles u bodu L2: rotující soustavu Země–Slunce, pseudopotenciál zahrnující centrifugální a Coriolisovy členy, linearizovanou analýzu stability u bodu L2 a pohybové rovnice dvou provázaných tuhých těles obíhajících v této soustavě. Referenční zdroje: Szebehely 1967, Murray & Dermott 1999. Mimo rozsah v0.1.

## §2. Závazky kanonické úrovně

Tato dokumentační úroveň se drží standardu **lege artis**. Stejná disciplína, která se uplatňuje v projektu `lege-artis/fourier` (viz `fourier/docs/canonical/en/00-introduction.md` §2), platí i zde:

- Každá pohybová rovnice je v této úrovni **odvozena z prvních principů**. Neříkáme pouze »moment síly z gravitačního gradientu je $3 G m / r^3 \cdot \hat{r} \times (\mathbf{I} \hat{r})$, protože to říká Hughes 2004« — ukazujeme multipolový rozvoj vzájemného gravitačního potenciálu, identifikujeme slapový člen, odvozujeme příslušný moment sil ze vztahu $\boldsymbol{\tau} = -\nabla_\theta U$ a ověřujeme, že se limit správně redukuje.
- Každé tvrzení odkazuje na primární zdroj. Viz `shared/reference-bibliography/refs.bib`. Pro formalismus analytické mechaniky: Thorne & Blandford 2017, Goldstein 2002, Arnold 1989, Marsden & Ratiu 1999, Landau–Lifšic, svazek I. Pro doménově specifickou dynamiku postoje kosmické lodě: Hughes 2004. Pro gravitační potenciály rozšířených těles: Murray & Dermott 1999, Jackson 1999 (multipolový formalismus přenesený z elektrodynamiky).
- **Wikipedie není primárním zdrojem.** Může sloužit jako rychlá kontrola, nikdy jako autoritativní pramen.
- **Numerické knihovny (NumPy, SciPy atd.) jsou orákula, nikoli zdrojový materiál.** Porovnání s `scipy.integrate.solve_ivp` slouží jako křížová validace vůči nezávislé implementaci; kanonické rovnice zůstávají pramenem pravdy.
- Matematický obsah existuje ve **dvou paralelních formách**: próza s rovnicemi v těchto souborech `docs/canonical/cs/NN-*.md` a strojově kontrolovatelné soubory rovnic v `shared/canonical-equations/*.{md,tex}` (připravováno pro kapitoly 02–05 v0.1).
- Dokumentace inženýrské úrovně (`docs/engineer/en/`, připravovaná) překládá tento obsah pro praktikanty, ale nikdy s ním není v rozporu. Pokud budoucí tvrzení inženýrské úrovně odporuje této kanonické úrovni, kanonická úroveň má přednost; prosíme o nahlášení jako issue.

## §3. Proč kanonická úroveň vůbec?

První verze na inženýrské úrovni v adresáři `experiments/jwst-l2-first-cut/` je postavena na Newtonových–Eulerových rovnicích vyjádřených přímo v inerciálních kartézských souřadnicích:

$$
\dot{\mathbf{p}}_i = \mathbf{F}_i, \qquad \dot{\mathbf{L}}_i = \boldsymbol{\tau}_i \quad \text{pro } i = A, B
$$

kde je úhlová rychlost přenášena v tělesové souřadnicové soustavě pomocí tenzoru setrvačnosti:

$$
\mathbf{I}\, \dot{\boldsymbol{\omega}} + \boldsymbol{\omega} \times (\mathbf{I}\, \boldsymbol{\omega}) = \boldsymbol{\tau}_{\text{body}}.
$$

To je správně a zákony zachování jsou počítány na strojovou přesnost (první verze dosahuje $|dE/E_0| = 6.2 \times 10^{-12}$ a $|dL/L_0| = 2.1 \times 10^{-10}$ při 600sekundové integraci s $dt = 0.05$ s). Nicméně v inženýrské úrovni není odvozeno — je tam *konstatováno* a validováno vůči analytickým limitům (volný (bezsilový) symetrický setrvačník atd.).

Čtenář kanonické úrovně se ptá: *proč* je pravá strana těchto rovnic taková, jaká je? Odkud pochází člen s vektorovým součinem úhlové rychlosti a tenzoru setrvačnosti? Proč je moment síly z gravitačního gradientu $3 G m / r^3 \cdot \hat{r} \times (\mathbf{I} \hat{r})$ v prvním řádu poměru rozměru tělesa k vzdálenosti, a co je člen vyššího řádu? Proč jsou Newtonovy–Eulerovy rovnice správné v inerciálních souřadnicích, ale Eulerovy rovnice správné v tělesové souřadnicové soustavě?

Kanonická úroveň odpovídá na tyto otázky postupně, počínaje variačním principem. Konfiguraci dvou tuhých těles pojímáme jako bod konfigurační variety $SE(3) \times SE(3)$ (speciální euklidovská grupa popisující tuhý pohyb, dvakrát), zapisujeme Lagrangeovu funkci $L = T - V$ na jejím tečném svazku, odvozujeme Eulerovy–Lagrangeovy rovnice a specializujeme. Newtonovy–Eulerovy rovnice inženýrské úrovně vycházejí jako konkrétní souřadnicový výraz tohoto; Eulerovy rovnice vycházejí z Euler-Poincarého redukce pomocí levé invariance kinetické energie vůči rotacím tělesa.

Výhoda tohoto přístupu není pouze pedagogická. Vynořují se tři konkrétní přínosy:

1. **Zákony zachování plynou ze struktury.** Zachování energie je invariance vůči časové translaci; zachování celkové hybnosti a momentu hybnosti jsou invariance vůči prostorové translaci, resp. rotaci. Noetherův teorém toto zpřístupňuje na úrovni Lagrangeovy funkce, ještě předtím, než se problému dotkne jakýkoli numerický integrátor.
2. **Redukované rovnice jsou jednodušší než ty získané redukcí z kartézských souřadnic.** Eulerovy dynamické rovnice v tělesové souřadnicové soustavě ($\mathbf{I}\,\dot{\boldsymbol{\omega}} + \boldsymbol{\omega} \times \mathbf{I}\,\boldsymbol{\omega} = \boldsymbol{\tau}$) jsou strukturálním důsledkem levé invariance kinetické energie na $SO(3)$, nikoli souřadnicovým trikem. Po tomto odvození se okamžitě zobecňují na jiné Lieovy grupy (např. $SE(3)$ pro těleso, které se zároveň pohybuje translačně, s propojenými rotační a translační dynamikou).
3. **Symplektické integrátory jsou nasnadě.** Inženýrská úroveň používá RK4, který je nesymplektický. Variační formulace kanonické úrovně navrhuje variační integrátory na Lieových grupách (Hairer–Lubich–Wanner 2006) jako strukturu zachovávající upgrade. V v0.1 je neimplementujeme; objasňujeme, co by dělaly a proč.

## §4. Požadavek odvozování z prvních principů — v konkrétních pojmech

Pro každou kapitolu v této úrovni:

- Každá rovnice má buď odvození v samotné kapitole, nebo citaci primárního zdroje, kde se odvození nachází. Žádné rovnice »konstatované bez důkazu«.
- Každé tvrzení o fyzikální interpretaci je ukotveno: »tento člen je centrifugální pseudosíla pociťovaná v rotující tělesové souřadnicové soustavě« musí být buď odvozeno, nebo citováno.
- Numerické příklady (vypočtené hodnoty, např. anizotropie tenzoru setrvačnosti tělesa podobného JWST a tělesa sondy) se odkazují zpět na kód inženýrské úrovně, který je produkuje, s citovaným rozsahem řádků.
- Validace vůči analytickým limitům je požadavkem každé kapitoly (Kapitola 06): každá kanonická rovnice musí být testována na alespoň jedné konfiguraci s uzavřeným analytickým řešením a nahlášen rozdíl.

Tato disciplína odráží akceptační kritéria C1–C4 projektu `lege-artis/fourier` (viz `_config/LEGE-ARTIS-LANGUAGE-DOCTRINE-v0.1.md` §4.4), specializovaná pro doménu analytické mechaniky. Kritéria jsou:

| Kritérium | Interpretace pro JWST-L2 |
|---|---|
| **C1** Pokrytí testy odpovídající učebnici | Kapitola 06 testuje každou kanonickou rovnici na ≥3 konfiguracích analytického limitu (volný (bezsilový) symetrický setrvačník; binární soustava bodových hmot — Kepler; volná precese v centrálním poli). **Plný katalog:** `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md` — Třída A (přesná uzavřená řešení) a Třída B (konfigurace okrajových případů), každá s matematickým nastavením, referenčním ukotvením, očekávanou numerickou tolerancí a umístěním testovacího souboru. |
| **C2** Bitová shoda mezi implementacemi (kde je aplikovatelné) | Až v0.2 přinese alternativní implementaci (např. symplektický variační integrátor na Lieových grupách vedle stávajícího RK4), bude porovnání implementací s tolerancí ULP × √N požadavkem CI. |
| **C3** Citace učebnic v záhlavích zdrojového kódu | Kód inženýrské úrovně (`dynamics.py`, `geometries.py`) nese v záhlavích citace na sekce kanonické úrovně, auditované programem grep před sestavením. (Připravováno pro v0.2 inženýrské úrovně.) |
| **C4** Empirický audit složitosti $O(N \log N)$ | Není přímo aplikovatelné na integraci ODE (která je $O(N \cdot N_{\mathrm{steps}})$). Nahrazeno: empirický audit chyby zachování vůči kroku (očekáváno $O(dt^4)$ pro drift energie s RK4). Viz `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md` §3.6 pro nastavení testu. |

Disciplína kanonické úrovně rozlišuje dvě doplňující se třídy testů, přičemž obě jsou nutné pro poctivou validaci modelu:

- **Třída A — Přesná uzavřená analytická řešení** jako orákula. V libovolném čase $t$ je analytická hodnota pravdou; vzdálenost numerické hodnoty od ní je měřitelná chyba. Odpovídá na otázku *»řeší můj integrátor rovnice správně?«*
- **Třída B — Konfigurace okrajových případů** testující doménu platnosti modelu. Odpovídá na otázku *»řeším správné rovnice?«*

Bez třídy B může být model numericky přesný, avšak fyzikálně chybný. Bez třídy A platí opak. Plný katalog testovacích konfigurací a okrajových případů v0.1 se nachází v `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md`; každý záznam se rozrůstá v odvozovací text Kapitoly 06, až tato kapitola bude hotova.

## §5. Doporučená čtecí cesta

Pro čtenáře se silným fyzikálním základem a předchozí zkušeností s Lagrangeovým formalismem:

1. Tento úvod (jste zde)
2. Kapitola 01 — Konfigurační varieta a variační princip (následující soubor)
3. Kapitola 02 — Kinetická energie na $SE(3) \times SE(3)$
4. Kapitola 03 — Vzájemný gravitační potenciál a multipolový rozvoj
5. Kapitola 04 — Eulerovy–Lagrangeovy rovnice a Euler-Poincarého redukce
6. Kapitola 05 — Eulerovy rovnice v tělesové souřadnicové soustavě jako specializace
7. Kapitola 06 — Analytické limity a validační brány

Pro čtenáře bez předchozí zkušenosti s Lagrangeovým formalismem doporučujeme začít s Feynman, svazek II, kap. 19 (Princip nejmenší akce) a Goldstein, kap. 1–2, poté se vrátit sem. Thorne & Blandford, svazek I nabízí stejný výchozí materiál s vyšší hustotou.

Pro čtenáře, který se zajímá spíše o kód než o matematiku, viz inženýrská úroveň `docs/engineer/en/00-quick-start.md` (připravovaná) nebo existující vstupní body první verze `experiments/jwst-l2-first-cut/dynamics.py` a `run_first_example.py`. Shad-úroveň `SHAD-NARRATIVE.md` je narativní průvodce první verzí.

## §6. Co tato úroveň záměrně nepokrývá

Pro přehlednost a vymezení očekávání:

- **Obecně relativistické korekce.** Orbitální rychlost JWST $v \approx 0.3 \text{ km/s}$ vůči barycentru soustavy Země–Slunce; $v/c \approx 10^{-6}$; efekt rotační vazby v prvním relativistickém řádu je úměrný $(v/c)^2 \approx 10^{-12}$, což je hluboko pod jakoukoli současnou měřitelnou přesností rotačního stavu. Pracujeme v čisté Newtonově gravitaci.
- **Tlak slunečního záření na sluneční clonu JWST.** Jde o dominantní negravitační perturbaci skutečné mise JWST a zodpovídá za podstatnou část palivového rozpočtu na udržení dráhy kosmické lodě. Jde o reálnou perturbaci; v kanonické úrovni ji vynecháváme, protože (a) nejde o gravitační vazbu, (b) závisí na detailní optice povrchu, která není dostupná v kanonických rovnicích, a (c) první verze prototypu inženýrské úrovně se zaměřuje specificky na gravitační vazbu. Verze v0.2 této kanonické úrovně zkatalogizuje, co by se změnilo, kdyby perturbace byla zahrnuta.
- **Vnitřní vibrační módy obou těles.** Předpokládáme tuhá tělesa. Skutečný JWST má strukturální pružnost; její modelování vyžaduje vrstvu mechaniky kontinua mimo rámec tuhého tělesa. Mimo rozsah projektu lege-artis.
- **Magnetické a elektrostatické interakce.** JWST je s vysokou přesností nenabité a provozuje se v blízkém vakuu; tyto interakce jsou v prostředí bodu L2 zanedbatelné a jsou ignorovány bez dalšího komentáře.
- **Perturbace s dlouhou periodou od Měsíce, Venuše, Jupitera.** Relevance pro návrh mise je reálná; rozsah kanonické úrovně je omezen na bezprostřední gravitační vazbu dvou těles. Budoucí rozšíření v0.3 by je přidalo.

## §7. Křížové odkazy a doprovodné materiály

- **Doprovodná inženýrská úroveň:** `docs/engineer/en/` (připravovaná v0.2 inženýrských dokumentů).
- **Doprovodná Shad-úroveň:** `SHAD-NARRATIVE.md` v kořenovém adresáři projektu (průvodce zásilkou k narozeninám psaný inženýrskou narací).
- **Referenční bibliografie:** `shared/reference-bibliography/refs.bib`.
- **Soubory rovnic** (připravované): `shared/canonical-equations/lagrangian.md`, `…/euler-equations.md`, `…/gravity-gradient-torque.md`, `…/multipole-expansion.md`. Každý soubor rovnic bude zkříženě odkazován ze sekce kapitoly, která ho odvozuje.
- **Kód první verze, jejž tato úroveň odvozuje:** `geometries.py`, `dynamics.py`, `run_first_example.py`, `testcases.py`, `tests/`.

---

**Následující kapitola:** [01 — Konfigurační varieta a variační princip](01-konfiguracni-varieta.md)

**Konec 00-uvod.md**
