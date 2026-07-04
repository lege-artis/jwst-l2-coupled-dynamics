# Kapitola 2 — Kinetická energie na $SE(3) \times SE(3)$ v intrinzické formě

> **Předpokládané čtení.** Kapitola 01 (konfigurační varieta a variační princip). Znalost tenzoru setrvačnosti jako multilineárního objektu je užitečná, avšak je zde pro úplnost vybudována od základu.
>
> **Stav.** v0.1.

> **Poznámka ke způsobu a předpokladům odvození.** Tato CS kanonická kapitola (v0.1) je odvozena z anglické kanonické kapitoly `docs/canonical/en/02-kinetic-energy.md` (v0.1, včetně doplnění 2. kola: §2.5.1, §2.6.1, §2.7.1). Terminologie je ukotvena v Podolský 2018 *Teoretická mechanika ve třech knihách* (Matfyzpress). EN verze zůstává autoritativním pramenem při budoucích aktualizacích.

---

## §2.1 Co tato kapitola dělá

Kapitola 1 stanovila konfigurační varietu $Q = SE(3) \times SE(3)$ a identifikovala tečný svazek $TQ$, na němž žije Lagrangián. Tato kapitola explicitně zapisuje **kinetickou energii**, polovinu onoho Lagrangiánu, odvozenou z prvních principů integrací přes rozložení hmoty každého tuhého tělesa. Výsledkem je hladká funkce $T: TQ \to \mathbb{R}$, jež je kvadratická v rychlostech a, pohlíženo geometricky, definuje **Riemannovu metriku** na $TQ$. Polovina s potenciální energií (kapitola 03) konstrukci uzavře; sestavení plného Lagrangiánu a odvození pohybových rovnic následuje v kapitole 04.

Odvození je **intrinzické**: kinetická energie je fyzikální skalár, týž v každé soustavě. Vyjadřujeme ji v souřadnicích jen tehdy, kdy to výpočet vyžaduje (tělesová soustava se ukazuje jako nejjednodušší volba souřadnic pro rotační část, inerciální soustava pro translační). Toto rozlišení intrinzické / souřadnicové vyjádření je důležité, neboť přežívá všechny následující redukce a specializace.

## §2.2 Kinetická energie jednoho tuhého tělesa — z rozložení hmoty

Uvažme jediné tuhé těleso $A$ s rozložením hmoty $\rho_A(\mathbf{x})$ neseným na jisté kompaktní oblasti $\mathcal{B}_A \subset \mathbb{R}^3$. Konfigurace tělesa v inerciální soustavě je popsána (dle kapitoly 1, §1.2) polohou těžiště $\mathbf{R}_A \in \mathbb{R}^3$ a orientací $q_A \in SO(3)$. Hmotný bod tělesa označený svou polohou $\mathbf{x} \in \mathcal{B}_A$ v tělesové soustavě (braný vzhledem k těžišti tělesa) sídlí v inerciální soustavě v:

$$
\mathbf{r}_A(\mathbf{x}, t) = \mathbf{R}_A(t) + \mathbf{R}_q^A(t)\, \mathbf{x},
$$

kde $\mathbf{R}_q^A$ je rotační matice příslušná $q_A$. Derivováním podle času je inerciální rychlost onoho hmotného bodu:

$$
\mathbf{v}_A(\mathbf{x}, t) = \dot{\mathbf{R}}_A + \dot{\mathbf{R}}_q^A\, \mathbf{x} = \dot{\mathbf{R}}_A + \hat{\boldsymbol{\omega}}_A\, \mathbf{R}_q^A\, \mathbf{x},
$$

s užitím vztahu $\dot{\mathbf{R}}_q^A = \hat{\boldsymbol{\omega}}_A\, \mathbf{R}_q^A$ z kapitoly 1, §1.4, kde $\boldsymbol{\omega}_A \in \mathbb{R}^3$ je úhlová rychlost tělesa $A$ vyjádřená v inerciální soustavě.

Kinetická energie tělesa $A$ je objemovým integrálem $\tfrac{1}{2} \rho |\mathbf{v}|^2$ přes těleso:

$$
T_A = \frac{1}{2} \int_{\mathcal{B}_A} \rho_A(\mathbf{x})\, |\mathbf{v}_A(\mathbf{x}, t)|^2\, dV.
$$

Dosazením výrazu pro rychlost a rozvinutím čtverce:

$$
|\mathbf{v}_A|^2 = |\dot{\mathbf{R}}_A|^2 + 2\, \dot{\mathbf{R}}_A \cdot \big(\hat{\boldsymbol{\omega}}_A\, \mathbf{R}_q^A\, \mathbf{x}\big) + |\hat{\boldsymbol{\omega}}_A\, \mathbf{R}_q^A\, \mathbf{x}|^2.
$$

Integrujme člen po členu. První člen dává $|\dot{\mathbf{R}}_A|^2 \int_{\mathcal{B}_A} \rho_A\, dV = m_A\, |\dot{\mathbf{R}}_A|^2$, kde $m_A$ je celková hmotnost. Druhý člen obsahuje integrál $\int \rho_A\, \mathbf{x}\, dV$, jenž **identicky mizí**, neboť $\mathbf{x}$ je měřeno od těžiště: z definice $\int \rho_A\, \mathbf{x}\, dV = \mathbf{0}$. Toto je algebraický obsah **věty o těžišti** (Königova věta).

Třetí člen vyžaduje opatrnost. S užitím antisymetrie $\hat{\boldsymbol{\omega}}$ a identity $|\mathbf{a} \times \mathbf{b}|^2 = |\mathbf{a}|^2 |\mathbf{b}|^2 - (\mathbf{a} \cdot \mathbf{b})^2$, s $\hat{\boldsymbol{\omega}}\, \mathbf{y} = \boldsymbol{\omega} \times \mathbf{y}$:

$$
|\hat{\boldsymbol{\omega}}_A\, \mathbf{R}_q^A\, \mathbf{x}|^2 = |\boldsymbol{\omega}_A \times (\mathbf{R}_q^A\, \mathbf{x})|^2 = |\boldsymbol{\omega}_A|^2 |\mathbf{R}_q^A\, \mathbf{x}|^2 - (\boldsymbol{\omega}_A \cdot \mathbf{R}_q^A\, \mathbf{x})^2.
$$

Jelikož $\mathbf{R}_q^A$ je ortogonální, $|\mathbf{R}_q^A\, \mathbf{x}| = |\mathbf{x}|$ a $\boldsymbol{\omega}_A \cdot \mathbf{R}_q^A\, \mathbf{x} = (\mathbf{R}_q^A)^T \boldsymbol{\omega}_A \cdot \mathbf{x} = \boldsymbol{\Omega}_A \cdot \mathbf{x}$, kde $\boldsymbol{\Omega}_A = (\mathbf{R}_q^A)^T \boldsymbol{\omega}_A$ je úhlová rychlost v **tělesové soustavě** (zavedená v kapitole 1, §1.4). Tudíž:

$$
|\hat{\boldsymbol{\omega}}_A\, \mathbf{R}_q^A\, \mathbf{x}|^2 = |\boldsymbol{\omega}_A|^2 |\mathbf{x}|^2 - (\boldsymbol{\Omega}_A \cdot \mathbf{x})^2.
$$

Provedením integrace se třetí člen stává:

$$
\int_{\mathcal{B}_A} \rho_A \left[ |\boldsymbol{\omega}_A|^2 |\mathbf{x}|^2 - (\boldsymbol{\Omega}_A \cdot \mathbf{x})^2 \right] dV.
$$

Toto je rozpoznatelně $\boldsymbol{\Omega}_A^T\, \mathbf{I}_A\, \boldsymbol{\Omega}_A$, kde $\mathbf{I}_A$ je **tenzor setrvačnosti tělesa $A$** v jeho tělesové soustavě (přesně definovaný v §2.4 níže); s užitím $|\boldsymbol{\omega}_A| = |\boldsymbol{\Omega}_A|$ (ortogonální transformace zachovávají normu) identita platí v obou soustavách.

## §2.3 Königův rozklad

Sebráním členů z §2.2 se kinetická energie tělesa $A$ rozkládá jako:

$$
\boxed{
T_A \;=\; \underbrace{\tfrac{1}{2}\, m_A\, |\dot{\mathbf{R}}_A|^2}_{\textstyle T_A^{\mathrm{trans}}}
\;+\; \underbrace{\tfrac{1}{2}\, \boldsymbol{\Omega}_A^T\, \mathbf{I}_A\, \boldsymbol{\Omega}_A}_{\textstyle T_A^{\mathrm{rot}}}.
}
$$

Toto je **Königova věta** pro kinetickou energii tuhého tělesa: celková kinetická energie se rovná kinetické energii pohybu těžiště plus rotační kinetické energii kolem těžiště. Smíšený člen zmizel díky volbě těžiště pro $\mathbf{x}$.

Dva důsledky hodné okamžitého zdůraznění:

1. Translační a rotační část se **rozpojují na úrovni kinetické energie**. Jakákoli vazba mezi translací a rotací musí pocházet z potenciální energie (kapitola 03) — konkrétně z členu gravitačního gradientu, jenž závisí na poloze i orientaci.

2. Rotační část je přirozeně vyjádřena v **tělesové soustavě** ($\boldsymbol{\Omega}_A$, $\mathbf{I}_A$), nikoli v inerciální. Tenzor setrvačnosti v tělesové soustavě je na čase nezávislý (těleso je tuhé; jeho rozložení hmoty se vůči vlastní soustavě nemění). Táž veličina v inerciální soustavě, $\mathbf{R}_q^A\, \mathbf{I}_A (\mathbf{R}_q^A)^T$, je na čase závislá skrze rotační matici. Tato asymetrie je tím, proč jsou tělesové souřadnice upřednostňovány pro rotační dynamiku — není to konvence, je to jediná volba, jež udržuje tenzor setrvačnosti statický.

> **Hlavní odkaz.** Pro odvození Königovy věty je hlavním CS-jazykovým ukotvením Podolský 2018 §7.1, kde je kinetická energie tuhého tělesa rozložena na translační a rotační část pomocí téhož argumentu mizejícího smíšeného členu při volbě těžiště. **Doplňující anglojazyčný odkaz.** Landau-Lifšic, sv. I §32 odvozuje Königovu větu v těsně paralelní formě (jejich značení: $\mathbf{V}$ pro rychlost těžiště, $\boldsymbol{\Omega}$ pro tělesovou úhlovou rychlost; týž obsah).

## §2.4 Tenzor setrvačnosti

**Tenzor setrvačnosti** $\mathbf{I}_A$ tělesa $A$ v jeho tělesové soustavě je symetrická pozitivně definitní matice $3 \times 3$:

$$
(\mathbf{I}_A)_{ij} \;=\; \int_{\mathcal{B}_A} \rho_A(\mathbf{x}) \left(|\mathbf{x}|^2\, \delta_{ij} - x_i\, x_j\right) dV.
$$

Diagonální prvky jsou momenty setrvačnosti kolem každé tělesové osy; mimodiagonální prvky jsou deviační momenty. Podle spektrální věty pro symetrické matice připouští $\mathbf{I}_A$ ortonormální vlastní bázi: **hlavní osy** tělesa. V tělesové soustavě hlavních os je $\mathbf{I}_A = \mathrm{diag}(I_{xx,A}, I_{yy,A}, I_{zz,A})$ diagonální s nezápornými prvky.

**Steinerova věta.** Je-li těleso složeno z podsoučástí $\{C_k\}$ s tenzory setrvačnosti $\mathbf{I}_{C_k}$ kolem jejich vlastních těžišť v tělesových polohách $\mathbf{c}_k$, je tenzor setrvačnosti složeného tělesa kolem jeho těžiště:

$$
\mathbf{I}_A \;=\; \sum_k \left[ \mathbf{I}_{C_k} + m_{C_k} \left(|\mathbf{c}_k|^2\, \mathbb{I}_3 - \mathbf{c}_k \mathbf{c}_k^T\right) \right].
$$

Druhý člen v závorce je Steinerův příspěvek: jde o tenzor setrvačnosti hmotného bodu $m_{C_k}$ v poloze $\mathbf{c}_k$ vzhledem k těžišti složeného tělesa.

> **Křížový odkaz na kód.** `geometries.py::parallel_axis(I_centroid, mass, r_offset)` implementuje přesně tento vzorec pomocí vnějšího součinu $\mathbf{r}\mathbf{r}^T$. Primitiva setrvačnosti složeného tělesa `disk_inertia`, `cylinder_inertia`, `thin_rod_inertia`, `solid_cone_inertia` poskytují hodnoty $\mathbf{I}_{C_k}$ pro standardní tvary podsoučástí. První verze `geometries.py` vystavuje tři referenční tělesa užívaná napříč pracovními příklady kanonické úrovně: (i) **`make_jwst_like()`** — zjednodušený 4složkový kompozit ve stylu JWST (sluneční clona + servisní modul + výložník + primární zrcadlo), celková hmotnost 2450 kg, produkuje $\mathbf{I} = \mathrm{diag}(2{,}33 \times 10^4,\ 2{,}33 \times 10^4,\ 1{,}54 \times 10^4)\ \text{kg m}^2$ (anizotropie max/min 1,52, **protáhlý**, neboť rozptyl složek podél z žene Steinerův příspěvek k $I_\perp$ nad lokální $I_\parallel$); užíván jako orákulum C2 křížového testu pro testy konstrukce kompozitu. (ii) **`make_probe_like()`** — malé 2složkové válcovo-talířové těleso, celková hmotnost 100 kg, $\mathbf{I} = \mathrm{diag}(29{,}06,\ 29{,}06,\ 5{,}00)\ \text{kg m}^2$ (anizotropie 5,81, protáhlý). (iii) **`make_oblate_reference_body()`** — syntetické jednoválcové těleso konstrukčně rozměrované tak, aby produkovalo přesně $\mathbf{I} = \mathrm{diag}(2540,\ 2540,\ 3870)\ \text{kg m}^2$ (anizotropie 1,52, **zploštělý**); užíváno jako referenční těleso pro pracovní příklad nestability gravitačního gradientu v §3.7.2 a příklad citlivosti d$\mathbf{I}$/d$t$ v §2.5.1 níže. Viz `_handoffs/fortran-phase-1/INVESTIGATION-OQ-FORT-6.md` pro historii, proč v tomto projektu koexistují dvě tělesa s anizotropií 1,52 (jedno protáhlé, jedno zploštělé).

Vlastní hodnoty $\mathbf{I}_A$ (hlavní momenty) jsou reálné a nezáporné. **Trojúhelníková nerovnost momentů** $I_{xx} + I_{yy} \ge I_{zz}$ (a cyklicky) je strukturálním faktem o rozloženích hmoty v 3D prostoru; jakýkoli „tenzor setrvačnosti“, jenž ji porušuje, neodpovídá žádnému fyzikálnímu rozložení hmoty.

> **Hlavní odkaz.** Pro algebru tenzoru setrvačnosti je hlavním CS-jazykovým ukotvením Podolský 2018 §7.1 Rovn. (7.7)–(7.9) (definice tenzoru setrvačnosti, $L_i = I_{ij}\Omega_j$ a $T = \tfrac{1}{2} I_{ij}\Omega_i\Omega_j$), (7.11) směrový moment a (7.12) Steinerova věta. **Doplňující anglojazyčné odkazy.** Goldstein 2002 §5.3 (tenzor setrvačnosti jako kartézský tenzor 2. řádu); Landau-Lifšic, sv. I §32 (hlavní osy a momenty); Hughes 2004 Kap. 4 (Steinerova věta pro kosmické lodě, včetně vzorů výpočtu složeného tělesa v podstatě totožných s přístupem `geometries.py`). Thorne & Blandford 2017 §1.1.1 („Geometrický pohled na zákony fyziky“) a §1.3 („Tenzorová algebra bez souřadnicové soustavy“) zarámovávají tenzor setrvačnosti v jazyce nezávislém na souřadnicích jako vektorově-hodnotovou lineární funkci $\mathbf{I}(\,\cdot\,, \boldsymbol{\omega}) = \mathbf{J}$ (moment hybnosti jako obraz úhlové rychlosti pod $\mathbf{I}$); jde o geometrický doplněk k souřadnicově-složkovému odvození výše. T&B v tomto svazku *neodvozuje* tenzor setrvačnosti z integrálu rozložení hmoty — materiál mechaniky tuhého tělesa náleží doprovodnému svazku nezahrnutému v nahraném vydání.

## §2.5 Tělesová vs. inerciální úhlová rychlost: účetnictví

Vztah $\boldsymbol{\Omega} = \mathbf{R}_q^T \boldsymbol{\omega}$ mezi tělesovou a inerciální úhlovou rychlostí zavedený v kapitole 1, §1.4 má několik jemností hodných upřesnění před kapitolou 3.

V tělesové soustavě je tenzor setrvačnosti $\mathbf{I}_A$ **konstantní v čase** (těleso je tuhé). Rotační kinetická energie je:

$$
T_A^{\mathrm{rot}} = \tfrac{1}{2}\, \boldsymbol{\Omega}_A^T\, \mathbf{I}_A\, \boldsymbol{\Omega}_A.
$$

V inerciální soustavě je tenzor setrvačnosti v čase $t$ roven $\mathbf{I}_A^{\mathrm{inertial}}(t) = \mathbf{R}_q^A(t)\, \mathbf{I}_A\, (\mathbf{R}_q^A(t))^T$, jenž je **na čase závislý skrze $\mathbf{R}_q^A(t)$**. Rotační kinetická energie je pak:

$$
T_A^{\mathrm{rot}} = \tfrac{1}{2}\, \boldsymbol{\omega}_A^T\, \mathbf{I}_A^{\mathrm{inertial}}(t)\, \boldsymbol{\omega}_A.
$$

Algebraicky jde o totéž číslo (dosaďte $\boldsymbol{\omega}_A = \mathbf{R}_q^A\, \boldsymbol{\Omega}_A$ a časová závislost $\mathbf{I}_A^{\mathrm{inertial}}$ se vyruší proti rotaci $\boldsymbol{\omega}_A$).

**Moment hybnosti** má analogické účetnictví:

$$
\mathbf{L}_A^{\mathrm{body}} = \mathbf{I}_A\, \boldsymbol{\Omega}_A, \qquad
\mathbf{L}_A^{\mathrm{inertial}} = \mathbf{R}_q^A\, \mathbf{L}_A^{\mathrm{body}} = \mathbf{I}_A^{\mathrm{inertial}}(t)\, \boldsymbol{\omega}_A.
$$

**Inerciální** moment hybnosti je zachovávající se veličinou, jsou-li vnější momenty sil nulové. Tělesový moment hybnosti se obecně nezachovává; splňuje Eulerovy rovnice (odvozené v kapitole 5). Tato asymetrie zachování je přímým důsledkem toho, která soustava je inerciální.

> **Křížový odkaz na kód.** `dynamics.py::inertia_inertial_frame(I_body, R_q)` rotuje tělesový tenzor setrvačnosti do inerciální soustavy podobnostní transformací $\mathbf{R}\, \mathbf{I}\, \mathbf{R}^T$. `dynamics.py::total_angular_momentum` počítá inerciální $\mathbf{L}_A$ užité pro diagnostiku zachování, pak sčítá přes tělesa. První verze dosahuje $|d\mathbf{L}/\mathbf{L}_0| = 2{,}1 \times 10^{-10}$ za 600 s, což je zaokrouhlovací podlaha pro aritmetiku dvojnásobné přesnosti při tomto integračním horizontu.

### §2.5.1 Časový vývoj $\mathbf{I}_A^{\mathrm{inertial}}(t)$: mimodiagonální citlivost (most na inženýrskou úroveň)

Tělesový tenzor setrvačnosti $\mathbf{I}_A$ je konstantní; inerciální tenzor setrvačnosti $\mathbf{I}_A^{\mathrm{inertial}}(t) = \mathbf{R}_q^A(t)\, \mathbf{I}_A\, (\mathbf{R}_q^A(t))^T$ je na čase závislý skrze orientaci. Landau-Lifšic a učebnicová tradice to považují za zřejmé a jdou dál. Pro implementaci reálných výpočtů — reprezentaci výstupu, analýzu citlivosti, volbu velikosti integračního kroku — jsou explicitní **časová derivace** $d\mathbf{I}_A^{\mathrm{inertial}}/dt$ a chování jejích **mimodiagonálních prvků** přirozeným mostem od kanonického formalismu k inženýrské praxi. Následující výklad je proto inženýrským doplňkem k §2.5, vyžádaným direktivou projektového vlastníka 2026-05-23.

Derivováním $\mathbf{I}_A^{\mathrm{inertial}}(t)$ podle času a užitím $\dot{\mathbf{R}}_q^A = \widehat{\boldsymbol{\omega}_A^{\mathrm{inertial}}}\, \mathbf{R}_q^A$ z kapitoly 1, §1.4 (s $\boldsymbol{\omega}_A^{\mathrm{inertial}} = \mathbf{R}_q^A \boldsymbol{\Omega}_A$ inerciální úhlovou rychlostí tělesa $A$):

$$
\frac{d\mathbf{I}_A^{\mathrm{inertial}}}{dt}
\;=\;
\dot{\mathbf{R}}_q^A\, \mathbf{I}_A\, (\mathbf{R}_q^A)^T
\;+\;
\mathbf{R}_q^A\, \mathbf{I}_A\, (\dot{\mathbf{R}}_q^A)^T
\;=\;
\widehat{\boldsymbol{\omega}_A^{\mathrm{inertial}}}\, \mathbf{I}_A^{\mathrm{inertial}}
\;-\;
\mathbf{I}_A^{\mathrm{inertial}}\, \widehat{\boldsymbol{\omega}_A^{\mathrm{inertial}}},
$$

tj. maticový komutátor:

$$
\boxed{
\frac{d\mathbf{I}_A^{\mathrm{inertial}}}{dt}
\;=\;
\big[\,\widehat{\boldsymbol{\omega}_A^{\mathrm{inertial}}}\, ,\, \mathbf{I}_A^{\mathrm{inertial}}\,\big].
}
$$

Toto je **rotační transportní rovnice** pro symetrický tenzor $\mathbf{I}_A^{\mathrm{inertial}}$. Je analogem, pro tenzor, elementárního $d\mathbf{v}/dt = \boldsymbol{\omega} \times \mathbf{v}$ pro vektor pevný v tělese. Stopa $\mathbf{I}_A^{\mathrm{inertial}}$ je zachována (komutátory mají nulovou stopu), stejně jako spektrum vlastních hodnot (podobnostní transformace zachovávají vlastní hodnoty) — konzistentně s tím, že hlavní momenty jsou intrinzickými vlastnostmi tělesa.

**Složkový tvar, soustava hlavních os v okamžiku vyhodnocení.** Uvažme okamžik, v němž $\mathbf{R}_q^A = \mathbb{I}$, takže $\mathbf{I}_A^{\mathrm{inertial}} = \mathbf{I}_A = \mathrm{diag}(I_1, I_2, I_3)$ je diagonální v souřadnicích inerciální soustavy. S tělesovou úhlovou rychlostí $\boldsymbol{\omega} = (\omega_1, \omega_2, \omega_3)$ (rovnou inerciální úhlové rychlosti v tomto okamžiku) přímý výpočet $[\widehat{\boldsymbol{\omega}}, \mathbf{I}_A]_{ij}$ dává rychlosti mimodiagonálních prvků jako:

$$
\left.\frac{d (\mathbf{I}_A^{\mathrm{inertial}})_{ij}}{dt}\right|_{\mathbf{R}_q^A = \mathbb{I}}
\;=\;
\varepsilon_{ijk}\, \omega_k\, (I_i - I_j), \qquad i \neq j.
$$

Rychlosti diagonálních prvků v tomto okamžiku mizí: $(d\mathbf{I}/dt)_{ii} = 0$ (což musí platit, jsou-li hlavní momenty zachovány). Vzorec je antisymetrický v $(i,j)$, konzistentně se symetrií $\mathbf{I}_A^{\mathrm{inertial}}$ v každém čase.

**Citlivost, podle třídy symetrie tělesa.** Vzorec pro mimodiagonální rychlost okamžitě odhaluje, jak se symetrie tělesa propagují do časového vývoje inerciálního tenzoru setrvačnosti:

| Třída symetrie tělesa | Hlavní momenty | Mimodiagonální rychlosti při $\mathbf{R}_q^A = \mathbb{I}$ |
|---|---|---|
| Sférická | $I_1 = I_2 = I_3 \equiv I_0$ | Vše nulové. $\mathbf{I}^{\mathrm{inertial}} = I_0 \mathbb{I}_3$ je invariantní vůči rotaci; žádná citlivost na orientaci v jakékoli veličině, jež se váže na $\mathbf{I}^{\mathrm{inertial}}$. |
| Axisymetrická ($I_1 = I_2 \neq I_3$) | Jeden odlišný moment | Rychlost $(1,2)$ mizí; rychlosti $(1,3)$ a $(2,3)$ jsou nenulové jen tehdy, má-li $\boldsymbol{\omega}$ kolmé složky $\omega_2$ resp. $\omega_1$. |
| Asymetrická (všechny odlišné) | Tři odlišné momenty | Všech šest mimodiagonálních rychlostí může být současně nenulových. |

**Numerická kontrola zdravého rozumu, zploštělé referenční těleso.** Užijme syntetické zploštělé axisymetrické referenční těleso `make_oblate_reference_body()` (viz křížový odkaz na kód v kapitole 2 §2.4 + řešení OQ-FORT-6 2026-05-24): $\mathbf{I}_A = \mathrm{diag}(2540, 2540, 3870)\ \text{kg m}^2$ (zploštělé, s $I_\parallel = I_3 = 3870 > I_\perp = I_1 = I_2 = 2540$, takže $I_\parallel / I_\perp = 1{,}52$). Pro pomalu rotující kosmickou loď s $\boldsymbol{\omega} = (0{,}01, 0, 0)\ \text{rad/s}$ (perioda rotace $T_{\mathrm{rot}} = 2\pi/|\boldsymbol{\omega}| = 628\ \text{s}$) výše uvedený vzorec předpovídá:
- $(d\mathbf{I}^{\mathrm{inertial}}/dt)_{2,3} = \omega_1 (I_3 - I_2) = 0{,}01 \cdot 1330 = 13{,}3\ \text{kg m}^2 / \text{s}$,

s ostatními mimodiagonálními rychlostmi nulovými díky axisymetrické struktuře. Čas, za nějž tento mimodiagonální prvek dosáhne typické diagonální velikosti $\sim 2540\ \text{kg m}^2$, je $\tau_{\mathrm{off}} \sim 2540 / 13{,}3 \approx 190\ \text{s} \approx T_{\mathrm{rot}}/3{,}3$ — v rámci integračního okna první verze $600\ \text{s}$ mimodiagonální prvek vykreslí přibližně jeden plný oscilační cyklus, amplitudou srovnatelný s diagonálními prvky. Mimodiagonální prvky **jsou** orientací tělesa, reprezentovanou ve složkách inerciální soustavy.

> **Dvě referenční tělesa, záměrně.** Práce kanonické úrovně užívá dvě komplementární axisymetrická referenční tělesa. (i) `make_oblate_reference_body()` — syntetické jednoválcové těleso s $\mathbf{I} = \mathrm{diag}(2540, 2540, 3870)$ konstrukčně, **zploštělé** ($I_\parallel > I_\perp$); užité zde a v kapitole 03 §3.7.2 pro numeriku pracovních příkladů. (ii) `make_jwst_like()` — zjednodušený 4složkový kompozit ve tvaru JWST (sluneční clona + servisní modul + výložník + primární zrcadlo) s $\mathbf{I} = \mathrm{diag}(2{,}33 \times 10^4, 2{,}33 \times 10^4, 1{,}54 \times 10^4)$, **protáhlý** ($I_\parallel < I_\perp$) díky Steinerovu příspěvku z rozptylu složek podél z; užitý pro demonstrace konstrukce kompozitu a jako orákulum C2 křížového testu. Kvalitativní srovnání obou je přirozenou křížovou validací: v librační rovnováze gravitačního gradientu je zploštělé těleso nestabilní (§3.7.2) a protáhlé těleso stabilní (§3.8.3). Viz `_handoffs/fortran-phase-1/INVESTIGATION-OQ-FORT-6.md` pro historii.

**Most k citlivosti momentu síly z gravitačního gradientu.** Tělesový moment síly z gravitačního gradientu odvozený v kapitole 03 §3.7,
$\boldsymbol{\tau}_A^{\mathrm{body}} = (3 G m_B / |\boldsymbol{\rho}|^3)\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \times (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A})$, závisí pouze na tělesových veličinách (konstantním $\mathbf{I}_A$ a orientačně závislém $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$). Inerciální pohled, ekvivalentní a užitečný pro reprezentaci výstupu, vyjadřuje týž moment síly jako
$\boldsymbol{\tau}_A^{\mathrm{inertial}} = (3 G m_B / |\boldsymbol{\rho}|^3)\, \hat{\boldsymbol{\rho}} \times (\mathbf{I}_A^{\mathrm{inertial}}(t) \cdot \hat{\boldsymbol{\rho}})$ s $\hat{\boldsymbol{\rho}}$ v inerciální soustavě. V této formě jsou časový vývoj momentu síly a jeho **citlivost na počáteční podmínky** přímo viditelné:
- **Citlivost na počáteční vzájemnou polohu** $\boldsymbol{\rho}(0)$: linearizovaně $\delta \boldsymbol{\tau}_A \propto \delta \boldsymbol{\rho} \cdot \nabla \boldsymbol{\tau}_A$, a Jacobiho matice $\nabla \boldsymbol{\tau}_A$ zahrnuje $\mathbf{I}_A^{\mathrm{inertial}}(0)$ — jejíž mimodiagonály kódují orientaci tělesa vůči inerciálním osám.
- **Citlivost na orientaci tělesa** $q_A(0)$: malé rotace $\delta q_A$ rotují $\mathbf{I}_A^{\mathrm{inertial}}$ o $\delta \mathbf{I}^{\mathrm{inertial}} = [\widehat{\delta q_A}, \mathbf{I}^{\mathrm{inertial}}]$ — táž komutátorová struktura jako $d\mathbf{I}/dt$. Odezva momentu síly je prvního řádu v této mimodiagonální perturbaci.
- **Citlivost sestavy dvou těles.** Celková dynamika se váže skrze časový vývoj orientace obou těles a jejich vzájemné $\boldsymbol{\rho}$. Pro dvě téměř sférická tělesa je dynamika přibližně rozpojená (moment síly z gravitačního gradientu mizí ve sférické limitě; §3.8.2 redukce). Pro dvě silně anizotropní tělesa v těsné vzájemné orientační sourodosti se efektivní anizotropie sestavy násobí a dynamika je odpovídajícím způsobem citlivější na počáteční podmínky — relevantním ukazatelem je **součin obou anizotropií** $(I_\parallel/I_\perp)_A \cdot (I_\parallel/I_\perp)_B$ vážený sourodostí $\hat{\boldsymbol{\rho}}$.

**Důsledek pro integrátor na inženýrské úrovni.** Mimodiagonální časové měřítko $\tau_{\mathrm{off}}$ odvozené výše stanovuje dolní mez rozlišeného časového měřítka rotační dynamiky. Aby RK4 tyto oscilace přesně sledoval, musí integrační krok $dt$ splňovat $dt \ll \tau_{\mathrm{off}}$. Konfigurace první verze užívá $dt = 0{,}05\ \text{s}$ proti $\tau_{\mathrm{off}} \approx 190\ \text{s}$, což dává zhruba čtyři řády rezervy — pohodlné. Pro těleso s extrémní anizotropií nebo velmi rychlou rotaci se rezerva zmenšuje a požadavek na velikost kroku se zpřísňuje; jde o druh praktického výpočetního omezení, jež most na inženýrskou úroveň zviditelňuje a jež samotný kanonický formalismus neukazuje.

> **Doplňující anglojazyčné odkazy.** Rotační transportní rovnice $d\mathbf{T}/dt = [\widehat{\boldsymbol{\omega}}, \mathbf{T}]$ pro symetrický tenzor v rotující soustavě je standardní; viz Hughes 2004 §4 (dynamika postoje kosmické lodě, formulovaná specificky pro tenzor setrvačnosti) a Goldstein 2002 §4.7-§4.10 (rychlosti změny v rotujících soustavách, obecný tenzorový případ). Analýza mimodiagonální citlivosti zde uvedená v pojmech tříd tělesové symetrie je implicitní v Hughes Kap. 8 (stabilizace gravitačním gradientem), avšak inženýrským publikem se přirozeněji čte, je-li zapsána explicitně.

> **Hlavní odkaz.** Podolský 2018 Kapitola 6 §6.3 odvozuje tutéž transportní rovnici v souřadnicově nezávislém tvaru, se zarámovaným výsledkem Rovn. (6.8): $\left(\frac{d\mathbf{w}}{dt}\right)\!\big|_{\mathrm{prostor}} = \left(\frac{d\mathbf{w}}{dt}\right)\!\big|_{\mathrm{t\check{e}leso}} + \boldsymbol{\Omega} \times \mathbf{w}$ pro libovolný vektor $\mathbf{w}$. Aplikováno sloupec po sloupci na sloupce $\mathbf{I}^{\mathrm{inertial}}$ (z nichž každý je vektorem pevným v tělese vyjádřeným v inerciálních souřadnicích) dává vektorový tvar (6.8) maticově-komutátorový tvar zarámovaný výše. Pomocná algebra v Podolský Rovn. (6.2) $\boldsymbol{\Omega}' = \frac{dA}{dt} A^t$ (tělesová) a (6.10) $\boldsymbol{\Omega} = A^t \frac{dA}{dt}$ (prostorová) poskytuje explicitní trivializaci, jež činí rozlišení tělesová-vs-inerciální soustava algebraicky transparentním.

## §2.6 Celková kinetická energie soustavy dvou těles

Kinetická energie dvou těles je součtem:

$$
T \;=\; T_A + T_B
\;=\; \tfrac{1}{2}\, m_A\, |\dot{\mathbf{R}}_A|^2 + \tfrac{1}{2}\, m_B\, |\dot{\mathbf{R}}_B|^2
+ \tfrac{1}{2}\, \boldsymbol{\Omega}_A^T\, \mathbf{I}_A\, \boldsymbol{\Omega}_A
+ \tfrac{1}{2}\, \boldsymbol{\Omega}_B^T\, \mathbf{I}_B\, \boldsymbol{\Omega}_B.
$$

Mezi tělesem $A$ a tělesem $B$ **nejsou žádné smíšené členy** v kinetické energii. Každé těleso přispívá svou vlastní translační + rotační kinetickou energií, nezávisle na druhém. Geometricky: kinetická energie na součinové varietě $TQ = T(SE(3)_A \times SE(3)_B) = TSE(3)_A \times TSE(3)_B$ se rozkládá jako součet kinetických energií na každém faktoru.

Toto **není zvláštní vlastnost** soustavy dvou těles; je to obecný fakt o kinetické energii mechanických soustav. Kinetická energie je součtem přes částice (nebo, ekvivalentně, přes podsystémy tuhých těles braných jako složené částice). Pojem „mechanické vazby“ mezi dvěma soustavami je v klasické mechanice vždy tvrzením o jejich potenciální energii, nikdy o kinetické.

(V Lagrangeově mechanice s **potenciály závislými na rychlosti** — např. nabitá částice v magnetickém poli, kde $L = \tfrac{1}{2} m \dot{\mathbf{x}}^2 + q \mathbf{A}(\mathbf{x}) \cdot \dot{\mathbf{x}}$ — člen lineární v rychlosti napodobuje kinetický smíšený člen, ale strukturálně jde o příspěvek potenciální energie. Pro naši čistě gravitační soustavu dvou těles žádný takový člen nevzniká a kinetická energie je čistým, v rychlosti kvadratickým součtem výše.)

Translační část můžeme přepsat pomocí **souřadnic těžiště + relativních** (konvence přijatá v kapitole 1, §1.3.1). Nechť $\mathbf{R}_{\mathrm{CoM}} = (m_A \mathbf{R}_A + m_B \mathbf{R}_B)/(m_A + m_B)$ a $\boldsymbol{\rho} = \mathbf{R}_B - \mathbf{R}_A$. Translační kinetická energie se stává:

$$
\tfrac{1}{2}\, m_A\, |\dot{\mathbf{R}}_A|^2 + \tfrac{1}{2}\, m_B\, |\dot{\mathbf{R}}_B|^2
= \tfrac{1}{2}\, M\, |\dot{\mathbf{R}}_{\mathrm{CoM}}|^2 + \tfrac{1}{2}\, \mu\, |\dot{\boldsymbol{\rho}}|^2,
$$

kde $M = m_A + m_B$ je celková hmotnost a $\mu = m_A m_B / M$ je **redukovaná hmotnost**. V těžišťové soustavě je $\dot{\mathbf{R}}_{\mathrm{CoM}} = \mathbf{0}$ identicky (toto je volba soustavy zamčená v kapitole 1) a translační kinetická energie se redukuje na člen relativního pohybu $\tfrac{1}{2}\, \mu\, |\dot{\boldsymbol{\rho}}|^2$.

V těžišťové soustavě je tedy celková kinetická energie:

$$
\boxed{
T_{\mathrm{CoM}} \;=\; \tfrac{1}{2}\, \mu\, |\dot{\boldsymbol{\rho}}|^2
\;+\; \tfrac{1}{2}\, \boldsymbol{\Omega}_A^T\, \mathbf{I}_A\, \boldsymbol{\Omega}_A
\;+\; \tfrac{1}{2}\, \boldsymbol{\Omega}_B^T\, \mathbf{I}_B\, \boldsymbol{\Omega}_B.
}
$$

Toto je tvar kinetické energie, jenž vstupuje do Lagrangiánu $L = T - V$ v kapitole 04. Změna proměnných $(\mathbf{R}_A, \mathbf{R}_B) \to (\mathbf{R}_{\mathrm{CoM}}, \boldsymbol{\rho})$ byla výše užita neformálně; následující podkapitola činí její status v kanonické mechanice explicitním.

### §2.6.1 Těžišťový rozklad jako symplektická bodová transformace

Výše uvedená změna proměnných je více než notačním pohodlím: je to **kanonická (symplektická) bodová transformace** na kotečném svazku konfigurační variety. Učinit toto explicitním je tím, co řeší OQ-2.1 ze seznamu OQ pracovní specifikace — upevňuje to formální status souřadnicové změny, již odvození kinetické energie užívá neformálně.

**Transformace konfiguračního prostoru.** Definujme lineární zobrazení na $\mathbb{R}^3 \times \mathbb{R}^3$:

$$
\begin{pmatrix} \mathbf{R}_{\mathrm{CoM}} \\ \boldsymbol{\rho} \end{pmatrix}
\;=\;
\mathbf{J}\, \begin{pmatrix} \mathbf{R}_A \\ \mathbf{R}_B \end{pmatrix},
\qquad
\mathbf{J} \;=\;
\begin{pmatrix} m_A/M & m_B/M \\ -\mathbb{I}_3 & +\mathbb{I}_3 \end{pmatrix}
\,\otimes\, \mathbb{I}_3
\;\equiv\;
\begin{pmatrix} m_A/M & m_B/M \\ -1 & +1 \end{pmatrix},
$$

kde druhá rovnost vypouští triviální faktor $\otimes \mathbb{I}_3$ a prvky matice jsou blokové skaláry působící nezávisle na každou prostorovou složku. Přímým výpočtem $\det \mathbf{J} = m_A/M + m_B/M = 1$. Transformace tedy zachovává objem na konfigurační varietě.

**Kotangenciální zdvih na $T^*Q$.** Standardní věta symplektické geometrie (Marsden & Ratiu 1999 §6.3; Arnold 1989 §3.16) říká, že jakýkoli difeomorfismus $\Phi: Q \to Q'$ se jednoznačně rozšiřuje na symplektický difeomorfismus $T^*\Phi: T^*Q \to T^*Q'$ — **kotangenciální zdvih** — pod nímž je kanonická 1-forma $\theta = p_i\, dq^i$ zachována. Pro lineární bodovou transformaci $Q' = \mathbf{J}\, Q$ působí zdvih na hybnosti inverzní transpozicí Jacobiho matice: $\mathbf{p}' = (\mathbf{J}^T)^{-1} \mathbf{p}$.

Výpočtem $\mathbf{J}^{-1}$ a její transpozice explicitně:

$$
\mathbf{J}^{-1} \;=\; \begin{pmatrix} 1 & -m_B/M \\ 1 & +m_A/M \end{pmatrix},
\qquad
(\mathbf{J}^T)^{-1} = (\mathbf{J}^{-1})^T \;=\; \begin{pmatrix} 1 & 1 \\ -m_B/M & m_A/M \end{pmatrix}.
$$

Indukovaná transformace hybnosti $(\mathbf{p}_A, \mathbf{p}_B) \to (\mathbf{P}_{\mathrm{CoM}}, \mathbf{p}_{\mathrm{rel}})$ je tedy:

$$
\boxed{
\mathbf{P}_{\mathrm{CoM}} \;=\; \mathbf{p}_A + \mathbf{p}_B, \qquad
\mathbf{p}_{\mathrm{rel}} \;=\; \frac{-m_B \mathbf{p}_A + m_A \mathbf{p}_B}{M}.
}
$$

**Ověření: $\mathbf{p}_{\mathrm{rel}} = \mu \dot{\boldsymbol{\rho}}$ podle očekávání.** Dosazením $\mathbf{p}_A = m_A \dot{\mathbf{R}}_A$, $\mathbf{p}_B = m_B \dot{\mathbf{R}}_B$:

$$
\mathbf{p}_{\mathrm{rel}}
\;=\; \frac{-m_B (m_A \dot{\mathbf{R}}_A) + m_A (m_B \dot{\mathbf{R}}_B)}{M}
\;=\; \frac{m_A m_B}{M}\, (\dot{\mathbf{R}}_B - \dot{\mathbf{R}}_A)
\;=\; \mu\, \dot{\boldsymbol{\rho}}.
$$

Hybnost konjugovaná k souřadnici relativní polohy je redukovaná hmotnost krát relativní rychlost, jak je požadováno.

**Kontrola symplektičnosti.** Výše uvedená věta o kotangenciálním zdvihu zaručuje, že $\theta$ je zachována, tudíž i $\omega = -d\theta$. Přímá kontrola to činí konkrétním: dosaďte inverzní transformace $\mathbf{p}_A = \mathbf{P}_{\mathrm{CoM}}\, m_A/M - \mathbf{p}_{\mathrm{rel}}$ a $\mathbf{p}_B = \mathbf{P}_{\mathrm{CoM}}\, m_B/M + \mathbf{p}_{\mathrm{rel}}$ do $d\mathbf{p}_A \wedge d\mathbf{R}_A + d\mathbf{p}_B \wedge d\mathbf{R}_B$, rozviňte pomocí $d\mathbf{R}_A = d\mathbf{R}_{\mathrm{CoM}} - (m_B/M)\, d\boldsymbol{\rho}$, $d\mathbf{R}_B = d\mathbf{R}_{\mathrm{CoM}} + (m_A/M)\, d\boldsymbol{\rho}$ a sebrání členů; smíšené členy se vyruší a výsledkem je $d\mathbf{P}_{\mathrm{CoM}} \wedge d\mathbf{R}_{\mathrm{CoM}} + d\mathbf{p}_{\mathrm{rel}} \wedge d\boldsymbol{\rho}$. 2-forma je zachována.

**Důsledek na straně Lagrangeova ($TQ$): kinetická metrika se blokově diagonalizuje.** Na tečném svazku užívá transformace rychlosti *tutéž* Jacobiho matici $\mathbf{J}$ jako transformace na straně konfigurace (lineární bodové transformace mají konstantní Jacobiho matice, takže strana rychlostí a konfigurace se transformují identicky). Metrika kinetické energie $\mathbf{g}_{\mathrm{old}} = \mathrm{diag}(m_A, m_B) \otimes \mathbb{I}_3$ se transformuje na:

$$
\mathbf{g}_{\mathrm{new}}
\;=\; (\mathbf{J}^{-1})^T\, \mathbf{g}_{\mathrm{old}}\, \mathbf{J}^{-1}
\;=\; \mathrm{diag}(M, \mu) \otimes \mathbb{I}_3.
$$

Přímým ověřením (se skalárním blokovým tvarem $\mathbf{J}^{-1}$ výše): off-bloky $(1,2)$ a $(2,1)$ mizí identicky a diagonální bloky vycházejí $M$ a $m_A m_B / M = \mu$. Toto je **algebraický obsah** struktury bez smíšených členů $T_{\mathrm{CoM}}$ v §2.6: Jacobiho matice těžišťového rozkladu je ortogonalizující transformací kinetické metriky vzhledem ke kanonickému skalárnímu součinu $(m_A, m_B)$. *Není* to náhoda; je to vynuceno geometrií pojmu těžiště.

**Perspektiva redukce.** Protože Lagrangián $L = T - V$ závisí na $\mathbf{R}_A, \mathbf{R}_B$ jen skrze jejich rozdíl $\boldsymbol{\rho}$ (translační invariance, kapitola 1 §1.6), je nová souřadnice $\mathbf{R}_{\mathrm{CoM}}$ cyklická: $\partial L / \partial \mathbf{R}_{\mathrm{CoM}} = 0$. Konjugovaná hybnost $\mathbf{P}_{\mathrm{CoM}}$ je tedy zachována (Noetherové věta aplikovaná na translační symetrii, jak bylo předjímáno v kapitole 1 v Tabulce symetrií). Těžišťovou souřadnici lze zcela eliminovat **symplektickou redukcí podle translační podgrupy** $\mathbb{R}^3 \subset SE(3)$ působící na $T^*(\mathbb{R}^3 \times \mathbb{R}^3)$, redukující translační fázový prostor z dimenze 12 na dimenzi 6 (relativní pohyb). Obecná teorie je v Marsden & Ratiu 1999 §1.2 + §11; pro naše účely fixace $\mathbf{P}_{\mathrm{CoM}} = \mathbf{0}$ a $\mathbf{R}_{\mathrm{CoM}} = \mathbf{0}$ (konvence klidové soustavy těžiště z kapitoly 1 §1.3.1) provádí redukci explicitně a ponechává jako aktivní translační stupně volnosti v zbylých odvozeních pouze proměnné relativního pohybu $(\boldsymbol{\rho}, \mathbf{p}_{\mathrm{rel}})$.

> **Doplňující anglojazyčné odkazy.** Marsden & Ratiu 1999 §6.3 (kotangenciální zdvih difeomorfismů — podkladová věta); §1.2 + §11 (symplektická redukce podle akce grupy — rámec, jenž z těžišťové symetrie činí počet dimenzí). Goldstein 2002 §9.4 pokrývá kanonické transformace v klasicko-mechanické tradici bez moderního geometrického zarámování; užitečné jako souřadnicově-složková křížová kontrola. Arnold 1989 §3.16 je stručným moderním pojednáním. **Hlavní odkaz.** Specifická bloková diagonalizace kinetické metriky dvou těles je v Landau-Lifšic, sv. I §13 (zarámování redukovanou hmotností, v souřadnicích); CS-jazyková tradice analytické mechaniky (Podolský) pojednává redukovanou hmotnost v rámci problému dvou těles týmž souřadnicovým rozkladem.

## §2.7 Kinetická energie jako Riemannova metrika na $TQ$ (geometrický pohled)

Až do tohoto bodu kapitola pracovala v souřadnicích. Intrinzické / souřadnicově nezávislé tvrzení toho, co jsme vybudovali, je následující.

**Riemannova metrika** na hladké varietě $M$ je hladké přiřazení, každému tečnému prostoru $T_p M$, pozitivně definitní symetrické bilineární formy $g_p: T_p M \times T_p M \to \mathbb{R}$. Kinetická energie mechanické soustavy, jejímž konfiguračním prostorem je $M$, definuje takovou metriku skrze $T = \tfrac{1}{2}\, g(\dot{q}, \dot{q})$. V souřadnicích se to stává známou kvadratickou formou $T = \tfrac{1}{2}\, g_{ij}(q)\, \dot{q}^i\, \dot{q}^j$.

Pro naši soustavu je $M = Q = SE(3)_A \times SE(3)_B$ a metrikou $g$ je ta definovaná kinetickou energií v §2.6. Struktura $g$ s redukovanou hmotností / setrvačností v těžišťové soustavě je:

| Tečný podprostor | Složka metriky |
|---|---|
| Rychlost relativní polohy $\dot{\boldsymbol{\rho}}$ | $\mu\, \mathbb{I}_3$ (euklidovská, škálovaná redukovanou hmotností) |
| Úhlová rychlost tělesa $A$ $\boldsymbol{\Omega}_A$ | $\mathbf{I}_A$ (tělesový tenzor setrvačnosti) |
| Úhlová rychlost tělesa $B$ $\boldsymbol{\Omega}_B$ | $\mathbf{I}_B$ (tělesový tenzor setrvačnosti) |
| Smíšené bloky | Vše nulové |

Metrika je **blokově diagonální**, neboť kinetická energie nemá smíšené členy (§2.6). Translační blok je **plochý** (euklidovský). Rotační bloky ploché nejsou — žijí na $SO(3)$, jež je Lieovou grupou s netriviální intrinzickou geometrií, a metrika na každém faktoru $SO(3)$ je **levo-invariantní**: hodnota metriky $g_q(\dot{q}, \dot{q})$ závisí pouze na tělesové úhlové rychlosti $\boldsymbol{\Omega}$, nikoli na aktuální orientaci $q$. Tato levá invariance je geometrickým obsahem tvrzení „tenzor setrvačnosti v tělesové soustavě je konstantní v čase“.

Levo-invariantní metriky na Lieových grupách jsou zdrojem struktury, již **Euler-Poincarého redukce** (kapitola 04) využívá. Lagrangián naší soustavy, má-li $T$ tento levo-invariantní tvar, se redukuje na pohybové rovnice na Lieově algebře $\mathfrak{so}(3) \oplus \mathfrak{so}(3) \oplus \mathbb{R}^3$ (tečný prostor relativní polohy) — šestirozměrném vektorovém prostoru — namísto přímo na 12rozměrné varietě $Q$. Tato redukce je tím, co produkuje **Eulerovy rovnice** pro rotaci tuhého tělesa v tělesové soustavě jako strukturální důsledek variačního principu, nikoli jako souřadnicový trik.

> **Doplňující anglojazyčné odkazy pro intrinzickou perspektivu.** Arnold 1989 §3.27 (tuhé těleso jako Lieova grupa); Marsden & Ratiu 1999 §9.1 (levo-invariantní metriky na Lieových grupách); Marsden & Ratiu §13 (Euler-Poincarého redukce, zde předjímaná, odvozená v kapitole 04). Pro čtenáře preferující elementárnější cestu Goldstein §5.3 pokrývá materiál tenzoru setrvačnosti / kinetické energie v souřadnicích a Landau-Lifšic, sv. I §32 ve stručné formě; oba dospívají k týmž souřadnicovým výrazům jako §2.6 bez geometrického zarámování.

### §2.7.1 Levá invariance rotační kinetické metriky na $SO(3)$: explicitní konstrukce

Hlavní text §2.7 tvrdil, že kinetická metrika na každém faktoru $SO(3)$ z $Q$ je **levo-invariantní**, a identifikoval to jako geometrický obsah časové nezávislosti tělesového tenzoru setrvačnosti. Tato podkapitola činí toto tvrzení konkrétním — co mechanicky znamená, že „metrika je levo-invariantní“ — a dokazuje je skrze explicitní identifikaci $g_R(\dot{q}_1, \dot{q}_2) = \langle \boldsymbol{\Omega}_1, \mathbf{I}\, \boldsymbol{\Omega}_2 \rangle_{\mathbb{R}^3}$ vyžádanou v OQ-2.3.

**Levá translace na $SO(3)$ a její diferenciál.** Pro každé $h \in SO(3)$ je **levá translace** o $h$ difeomorfismus
$$
L_h : SO(3) \to SO(3), \qquad L_h(R) \;=\; h R.
$$

Její diferenciál v $R$ je lineární zobrazení $(dL_h)_R : T_R SO(3) \to T_{hR} SO(3)$ mezi tečnými prostory. Pro tečný vektor $\dot{q} \in T_R SO(3)$ — reprezentovaný jako rychlost křivky procházející $R$ — je působení diferenciálu prostě maticovým násobením zleva:
$$
(dL_h)_R\, \dot{q} \;=\; h \dot{q},
$$
což je opět tečný vektor, nyní ukotvený v $h R$.

**Levo-invariantní vektorová pole.** Hladké vektorové pole $X$ na $SO(3)$ je **levo-invariantní**, jestliže $(dL_h)_R X(R) = X(hR)$ pro všechna $h, R \in SO(3)$. Ekvivalentně je $X$ určeno napříč $SO(3)$ svou hodnotou $X(\mathbb{I}_3) \in T_{\mathbb{I}_3} SO(3) = \mathfrak{so}(3)$ v identitě a propagováno do ostatních bodů levou translací. Prostor levo-invariantních vektorových polí je tedy izomorfní **Lieově algebře** $\mathfrak{so}(3)$ antisymetrických matic $3 \times 3$, jež je zase izomorfní $\mathbb{R}^3$ skrze operátor stříšky z kapitoly 1, §1.4. Oba izomorfismy jsou lineární a trojrozměrné; společně zdůvodňují zacházení s vektory úhlové rychlosti $\boldsymbol{\Omega} \in \mathbb{R}^3$ jako s přirozenými souřadnicemi na $\mathfrak{so}(3)$.

**Levá trivializace tečného vektoru.** Pro libovolné $\dot{q} \in T_R SO(3)$ je **levá trivializace** $\dot{q}$ prvkem $\widehat{\boldsymbol{\Omega}} \in \mathfrak{so}(3)$ získaným stažením $\dot{q}$ k identitě skrze $L_{R^{-1}}$:
$$
\widehat{\boldsymbol{\Omega}} \;\equiv\; (dL_{R^{-1}})_R\, \dot{q} \;=\; R^{-1} \dot{q} \;=\; R^T \dot{q}.
$$
Ekvivalentně $\dot{q} = R \widehat{\boldsymbol{\Omega}}$. Toto je přesně tělesový vztah z kapitoly 1, §1.4: tělesová úhlová rychlost $\boldsymbol{\Omega}$ je levou trivializací rychlosti orientace $\dot{q}$. Tělesová soustava je **levo-trivializovanou soustavou** na $SO(3)$.

**Kinetická metrika je levo-invariantní.** Definujme metriku na $SO(3)$ jako:
$$
\boxed{
g_R(\dot{q}_1, \dot{q}_2) \;\equiv\; \langle \boldsymbol{\Omega}_1, \mathbf{I}\, \boldsymbol{\Omega}_2 \rangle_{\mathbb{R}^3}, \qquad \widehat{\boldsymbol{\Omega}_i} = R^T \dot{q}_i,
}
$$
kde $\mathbf{I}$ je **konstantní** tělesový tenzor setrvačnosti (symetrická pozitivně definitní matice $3 \times 3$, nezávislá na $R$) a $\langle \cdot , \cdot \rangle_{\mathbb{R}^3}$ je standardní euklidovský skalární součin na $\mathbb{R}^3$.

K důkazu, že $g$ je levo-invariantní, vyhodnoťme v translovaném bodě $h R$ s translovanými tečnými vektory $(dL_h)_R \dot{q}_i = h \dot{q}_i$:
$$
g_{hR}\big((dL_h) \dot{q}_1,\, (dL_h) \dot{q}_2\big)
\;=\; \langle \boldsymbol{\Omega}'_1, \mathbf{I}\, \boldsymbol{\Omega}'_2 \rangle_{\mathbb{R}^3},
$$
kde $\widehat{\boldsymbol{\Omega}'_i} = (hR)^T (h \dot{q}_i) = R^T h^T h \dot{q}_i = R^T \dot{q}_i = \widehat{\boldsymbol{\Omega}_i}$, s užitím $h^T h = \mathbb{I}_3$ pro $h \in SO(3)$. Levé trivializace jsou levou translací nezměněny, tudíž:
$$
g_{hR}\big((dL_h) \dot{q}_1,\, (dL_h) \dot{q}_2\big)
\;=\; \langle \boldsymbol{\Omega}_1, \mathbf{I}\, \boldsymbol{\Omega}_2 \rangle_{\mathbb{R}^3}
\;=\; g_R(\dot{q}_1, \dot{q}_2). \qquad \blacksquare
$$

**Mechanická interpretace: tvrzení „$\mathbf{I}_{\mathrm{body}}$ je konstantní“ *je* tvrzení „$g$ je levo-invariantní“.** Nejde o dvě oddělené vlastnosti. To, že tělesový tenzor setrvačnosti $\mathbf{I}$ je v čase konstantní, znamená, že nezávisí na orientaci $R$; ekvivalentně, je-li metrika vyjádřena v tělesové ($\boldsymbol{\Omega}$) parametrizaci tečných vektorů, neobjeví se žádné na $R$ závislé koeficienty. Hodnota metriky $g_R(\dot{q}_1, \dot{q}_2) = \langle \boldsymbol{\Omega}_1, \mathbf{I}\, \boldsymbol{\Omega}_2 \rangle$ závisí na $R$ jen skrze trivializaci $\boldsymbol{\Omega}_i = R^T \dot{q}_i$ a tato závislost se vyruší pod současnou levou translací $R$ a tečných vektorů. To vyrušení je algebraickým obsahem levé invariance.

**Kontrast: pravá trivializace a inerciální popis.** Zrcadlová konstrukce užívá **pravou translaci** $R_h(R) = R h$ namísto levé. Jejím diferenciálem je násobení zprava: $(dR_h)_R \dot{q} = \dot{q} h$. **Pravá trivializace** $\dot{q}$ extrahuje $\widehat{\boldsymbol{\omega}} = \dot{q} R^{-1}$, **inerciální** úhlovou rychlost z kapitoly 1, §1.4. Táž metrika vyjádřená v této trivializaci:
$$
g_R(\dot{q}_1, \dot{q}_2) \;=\; \langle \boldsymbol{\omega}_1, \mathbf{I}_R^{\mathrm{inertial}}\, \boldsymbol{\omega}_2 \rangle_{\mathbb{R}^3}, \qquad \mathbf{I}_R^{\mathrm{inertial}} = R\, \mathbf{I}\, R^T,
$$
užívá na $R$ závislý **inerciální tenzor setrvačnosti** $\mathbf{I}_R^{\mathrm{inertial}}$ — přesně na čase závislý tenzor, jehož vývoj jsme studovali v §2.5.1. Táž metrika je nyní zjevně **pravo-invariantní** (pravo-translační analog levé invariance, s inerciální soustavou hrající roli, již výše hrála tělesová), avšak *není* levo-invariantní: změna $R$ mění matici $\mathbf{I}_R^{\mathrm{inertial}}$. Táž fyzikální metrika připouští dvě komplementární trivializace; *tělesová* soustava je **levo**-invariantní soustavou a *inerciální* soustava je **pravo**-invariantní soustavou.

Toto je geometrickým zdrojem asymetrie již pozorované v §2.5: tenzor setrvačnosti je konstantní v tělesové soustavě (levo-invariantní soustavě), avšak moment hybnosti se zachovává v inerciální soustavě (pravo-ekvivariantní). Konstantnost $\mathbf{I}_{\mathrm{body}}$ pochází z levé invariance metriky; zachování inerciálního $\mathbf{L}$ pochází z toho, že pravá akce (rotace inerciální soustavy) je Noetherovou symetrií Lagrangiánu soustavy.

**Proč na levé invarianci záleží: náhled na Euler-Poincarého redukci.** Lagrangián, jenž je levo-invariantní na Lieově grupě $G$ — tj. $L(g, \dot{g}) = L(\mathbb{I}, g^{-1} \dot{g})$ pro všechna $g$ — se redukuje na Lagrangián $\ell(\boldsymbol{\Omega})$ na Lieově algebře $\mathfrak{g} = T_{\mathbb{I}} G$, kde $\boldsymbol{\Omega} = g^{-1} \dot{g}$ je levou trivializací. Eulerovy–Lagrangeovy rovnice na $TG$ se pak transformují na **Euler-Poincarého rovnice** na $\mathfrak{g}^*$:
$$
\frac{d}{dt} \frac{\partial \ell}{\partial \boldsymbol{\Omega}} \;=\; \mathrm{ad}^*_{\boldsymbol{\Omega}}\, \frac{\partial \ell}{\partial \boldsymbol{\Omega}}.
$$
Pro rotační kinetickou energii $\ell(\boldsymbol{\Omega}) = \tfrac{1}{2} \boldsymbol{\Omega}^T \mathbf{I}\, \boldsymbol{\Omega}$ na $\mathfrak{so}(3) \cong \mathbb{R}^3$ se toto specializuje (s $\partial \ell / \partial \boldsymbol{\Omega} = \mathbf{I}\, \boldsymbol{\Omega}$ a $\mathrm{ad}^*_{\boldsymbol{\Omega}}\, v = -\boldsymbol{\Omega} \times v$ pro $\mathfrak{so}(3)$ ztotožněnou s $\mathbb{R}^3$) na:
$$
\frac{d}{dt}(\mathbf{I}\, \boldsymbol{\Omega}) \;=\; -\boldsymbol{\Omega} \times (\mathbf{I}\, \boldsymbol{\Omega}),
$$
což se přeskupuje na **Eulerovy rovnice** $\mathbf{I}\, \dot{\boldsymbol{\Omega}} + \boldsymbol{\Omega} \times (\mathbf{I}\, \boldsymbol{\Omega}) = \mathbf{0}$ (bez momentu sil; s vnějším momentem sil na pravé straně v obecném případě). Tělesové Eulerovy rovnice jsou tedy **strukturálním důsledkem** levé invariance kinetické metriky, nikoli souřadnicovým trikem. Kapitola 04 provádí tuto redukci podrobně; kapitola 05 se specializuje na tělesové Eulerovy rovnice explicitně.

**Poznámka k inženýrské úrovni.** Implementace první verze `dynamics.py` parametrizuje orientaci jednotkovým kvaternionem (mapou na $SO(3)$) a propaguje tělesovou úhlovou rychlost $\boldsymbol{\Omega}_A$ skrze Eulerovy rovnice (řádky 162–163), pak rekonstruuje orientaci kvaternionovou kinematikou (řádky 167–168). Toto je **levo-trivializovaná** integrace: dynamické proměnné $(\boldsymbol{\Omega}_A, q_A)$ jsou levou trivializací $\boldsymbol{\Omega}$ spárovanou s orientací $q$ a pohybové rovnice v těchto proměnných jsou časově autonomní (neboť $\mathbf{I}_A$ je konstantní). Alternativní pravo-trivializovaná integrace by propagovala inerciální $\boldsymbol{\omega}_A$ s na čase závislým $\mathbf{I}_A^{\mathrm{inertial}}(t)$ — ekvivalentní fyzikálním obsahem, avšak algebraicky nákladnější na krok, vyžadující aktualizaci orientačně závislé matice setrvačnosti při každém podkroku RK4. Levá trivializace je správnou volbou pro výpočetní efektivitu právě proto, že metrika je levo-invariantní.

> **Doplňující anglojazyčné odkazy.** Marsden & Ratiu 1999 §9.1 — kanonický výklad levo-invariantních metrik na Lieových grupách, s $SO(3)$ jako pracovním příkladem a Euler-Poincarého redukcí v §13 jako aplikací. Arnold 1989 §3.27 (tuhé těleso jako Lieova grupa) a Arnoldův dodatek o hamiltonovském formalismu pokrývají týž materiál. Pro těsnější moderní prezentaci Holm-Schmah-Stoica 2009 *Geometric Mechanics and Symmetry* §6 vede konstrukci od $SO(3)$ přes $SE(3)$ k obecným maticovým Lieovým grupám, užitečná příprava na redukci plného Lagrangiánu dvou těles v kapitole 04.

> **Hlavní odkaz.** Podolský 2018 Kapitola 6 §6.3 zavádí tenzor tělesové úhlové rychlosti skrze Rovn. (6.4) $\Omega'_{ij} \equiv \Omega(\mathbf{e}'_i, \mathbf{e}'_j)$ — složky tenzoru úhlové rychlosti vyhodnocené na bázových vektorech pevných v tělese (korotujících). Trivializační algebra v Rovn. (6.2) $\boldsymbol{\Omega}' = \frac{dA}{dt} A^t$ (tělesové složky) a (6.10) $\boldsymbol{\Omega} = A^t \frac{dA}{dt}$ (prostorové složky) je inženýrským vykreslením rozlišení levá-vs-pravá trivializace učiněného zde geometricky: Podolského „$\boldsymbol{\Omega}'$ v korotující bázi“ je levo-trivializovaná tělesová úhlová rychlost, jeho „$\boldsymbol{\Omega}$ v prostorové bázi“ je pravo-trivializovaná inerciální úhlová rychlost a ekvivalence obou fyzikálních popisů odpovídá téže metrice připouštějící levo- i pravo-invariantní vykreslení na různých stranách $SO(3)$. Podolský nevolá jazyk Lieových grup explicitně, avšak algebra, již odvozuje, je tělesová/prostorová trivializace v souřadnicích a je tím, co moderní Marsden-Ratiu pojednání abstrahuje do souřadnicově nezávislé formy.

## §2.8 Co je stanoveno na konci této kapitoly

Čtenář, jenž prošel §2.1–§2.7, má:

- Kinetickou energii $T$ soustavy dvou tuhých těles jako explicitní funkci na tečném svazku $TQ$, odvozenou z prvních principů integrací přes rozložení hmoty každého tělesa.
- Königův rozklad na translační + rotační část, s větou o těžišti jako algebraickým základem toho, proč nevzniká smíšený člen.
- Tenzor setrvačnosti $\mathbf{I}_A$ jako symetrickou pozitivně definitní matici $3 \times 3$ se Steinerovým pravidlem skládání, mapující přímo na konstrukci složeného tělesa `geometries.py` v kódu první verze inženýrské úrovně.
- Účetnictví tělesové vs. inerciální úhlové rychlosti, s explicitní identifikací toho, která soustava udržuje tenzor setrvačnosti statický (tělesová) vs. ve které soustavě se moment hybnosti zachovává při nepřítomnosti vnějšího momentu sil (inerciální).
- **Rotační transportní rovnici** $d\mathbf{I}^{\mathrm{inertial}}/dt = [\widehat{\boldsymbol{\omega}}, \mathbf{I}^{\mathrm{inertial}}]$ a citlivost jejích mimodiagonálních prvků na tělesovou symetrii a úhlovou rychlost (§2.5.1) — most na inženýrskou úroveň, jenž překládá kanonický formalismus do praktických tvrzení o citlivosti na počáteční podmínky, anizotropii tělesa a volbu velikosti integračního kroku.
- Tvar translační kinetické energie s redukovanou hmotností v těžišťové soustavě; celkové $T_{\mathrm{CoM}}$ jako vstupní bod do Lagrangiánu.
- **Těžišťový rozklad jako symplektickou bodovou transformaci** (§2.6.1) — změnu proměnných $(\mathbf{R}_A, \mathbf{R}_B) \to (\mathbf{R}_{\mathrm{CoM}}, \boldsymbol{\rho})$ s přesným statusem kanonické mechaniky jako kotangenciální zdvih lineární bodové transformace s jednotkovým determinantem Jacobiho matice, s indukovanou transformací hybnosti $\mathbf{P}_{\mathrm{CoM}} = \mathbf{p}_A + \mathbf{p}_B$, $\mathbf{p}_{\mathrm{rel}} = \mu \dot{\boldsymbol{\rho}}$ a blokovou diagonalizací kinetické metriky $\mathrm{diag}(M, \mu)$.
- Intrinzické geometrické tvrzení $T$ jako Riemannovy metriky na $TQ$, blokově diagonální, levo-invariantní na faktorech $SO(3)$. Toto je strukturální vlastnost, již Euler-Poincarého redukce (kapitola 04) využije.
- **Explicitní konstrukci** levé invariance pro rotační kinetickou metriku na každém faktoru $SO(3)$ (§2.7.1), se vzorcem $g_R(\dot{q}_1, \dot{q}_2) = \langle \boldsymbol{\Omega}_1, \mathbf{I}\, \boldsymbol{\Omega}_2 \rangle$ odvozeným a dokázaným levo-invariantním přímým výpočtem. Mechanická ekvivalence „$\mathbf{I}_{\mathrm{body}}$ konstantní $\Leftrightarrow$ $g$ levo-invariantní“ je identifikována; kontrast s inerciální (pravo-trivializovanou) parametrizací je učiněn explicitním; Euler-Poincarého redukce je předjímána jako důsledek levé invariance, jenž produkuje Eulerovy rovnice v kapitolách 04+05.

**Co dosud stanoveno NENÍ:**

- Potenciální energie $V$ explicitně. To vyžaduje vzájemnou gravitační interakci mezi rozšířenými tělesy, odvozenou multipolovým rozvojem. **Kapitola 03.**
- Plný Lagrangián $L = T - V$ a Eulerovy–Lagrangeovy pohybové rovnice na $TQ$. **Kapitola 04.**
- Tělesové Eulerovy rovnice $\mathbf{I}\dot{\boldsymbol{\Omega}} + \boldsymbol{\Omega} \times (\mathbf{I}\boldsymbol{\Omega}) = \boldsymbol{\tau}$ jako specializace Eulerových–Lagrangeových rovnic skrze Euler-Poincarého redukci. **Kapitola 05.**
- Validace vůči analytickým limitům — bezsilový symetrický setrvačník, nestabilita asymetrického setrvačníku, limit bodové hmoty Kepler atd. **Kapitola 06**, s plným katalogem zamčeným v `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md`.

---

## §2.9 Otevřené otázky v této kapitole

- **OQ-2.1.** ***Vyřešeno 2026-05-23*** — řešeno v nové podkapitole §2.6.1 („Těžišťový rozklad jako symplektická bodová transformace“). Změna proměnných $(\mathbf{R}_A, \mathbf{R}_B) \to (\mathbf{R}_{\mathrm{CoM}}, \boldsymbol{\rho})$ je lineární bodovou transformací s Jacobiho maticí $\mathbf{J}$ ($\det \mathbf{J} = 1$); její kotangenciální zdvih na $T^*(\mathbb{R}^3 \times \mathbb{R}^3)$ je symplektickou transformací mapující $(\mathbf{p}_A, \mathbf{p}_B)$ na $(\mathbf{P}_{\mathrm{CoM}} = \mathbf{p}_A + \mathbf{p}_B,\, \mathbf{p}_{\mathrm{rel}} = \mu \dot{\boldsymbol{\rho}})$; kanonická 2-forma je zachována (ověřeno přímo); metrika kinetické energie se blokově diagonalizuje na $\mathrm{diag}(M, \mu) \otimes \mathbb{I}_3$. Formální ukotvení: Marsden & Ratiu 1999 §6.3 (kotangenciální zdvih) + §1.2 + §11 (symplektická redukce). Na straně tečného svazku užívá souřadnicová změna rychlosti tutéž Jacobiho matici $\mathbf{J}$ a produkuje tutéž blokovou diagonalizaci kinetické metriky.
- **OQ-2.2.** ***Vyřešeno 2026-05-23*** — projektovým vlastníkem dodané T&B PDF zrevidováno (1552stránkové elektronické vydání). **Negativní nález pro Königovu větu (§2.3):** odvození Königova rozkladu se v T&B svazku *Modern Classical Physics* neobjevuje; téma mechaniky tuhého tělesa náleží doprovodnému svazku mechaniky nezahrnutému v nahraném vydání. Hlavním CS-jazykovým ukotvením §2.3 je Podolský §7.1; doplňujícím anglojazyčným zůstává Landau-Lifšic, sv. I §32. **Částečný nález pro tenzor setrvačnosti (§2.4):** T&B §1.1.1 + §1.3 (souřadnicově nezávislé geometrické zarámování $\mathbf{I}(\,\cdot\,, \boldsymbol{\omega}) = \mathbf{J}$) přidáno do referenčního bloku §2.4 jako geometrický doplněk k souřadnicovému odvození. Definice $\mathbf{I}$ integrálem rozložení hmoty v §2.4 si ponechává původní ukotvení na Goldstein/LL/Hughes.
- **OQ-2.3.** ***Vyřešeno 2026-05-23*** — řešeno v nové podkapitole §2.7.1 („Levá invariance rotační kinetické metriky na $SO(3)$: explicitní konstrukce“). Konkrétní výstup: levá translace $L_h$ a její diferenciál definovány explicitně; levo-invariantní vektorová pole ztotožněna s $\mathfrak{so}(3) \cong \mathbb{R}^3$; levá trivializace $\widehat{\boldsymbol{\Omega}} = R^T \dot{q}$ spárována s konvencí tělesové úhlové rychlosti z kapitoly 1 §1.4; metrika $g_R(\dot{q}_1, \dot{q}_2) = \langle \boldsymbol{\Omega}_1, \mathbf{I}\, \boldsymbol{\Omega}_2 \rangle$ (zarámovaná) ukázána levo-invariantní přímým algebraickým výpočtem; mechanická ekvivalence „$\mathbf{I}_{\mathrm{body}}$ konstantní $\Leftrightarrow$ $g$ levo-invariantní“ vyslovena a dokázána; kontrast s pravo-trivializovanou / inerciální formulací dán pro kontext; náhled na Euler-Poincarého redukci s explicitní specializací na Eulerovu rovnici; poznámka k inženýrské úrovni vážící levo-trivializovanou parametrizaci přímo na volbu integrace Eulerových rovnic v `dynamics.py`. Ukotvení: Marsden & Ratiu 1999 §9.1 + §13; Arnold 1989 §3.27.
- **OQ-2.4.** ***Vyřešeno 2026-05-23*** — most na inženýrskou úroveň pro formalismus tenzoru setrvačnosti zpracován v nové podkapitole §2.5.1 („Časový vývoj $\mathbf{I}_A^{\mathrm{inertial}}(t)$: mimodiagonální citlivost“). Zachycuje direktivu projektového vlastníka 2026-05-23: „Landau-Lifšic je zlatý standard, avšak pro implementaci reálných výpočtů by časová derivace tenzoru setrvačnosti — zejména jeho mimodiagonálních prvků — poskytla most k ilustraci toho, jak citlivá je dynamika na počáteční vzájemné polohy a jak osy a symetrie objektů a celé sestavy dvou těles vstupují.“ Výstupy: rotační transportní rovnice $d\mathbf{I}^{\mathrm{inertial}}/dt = [\widehat{\boldsymbol{\omega}}, \mathbf{I}^{\mathrm{inertial}}]$ s explicitním vzorcem pro mimodiagonální rychlost $(d\mathbf{I}/dt)_{ij} = \varepsilon_{ijk}\, \omega_k\, (I_i - I_j)$ v okamžiku $\mathbf{R}_q = \mathbb{I}$; tabulka citlivosti podle třídy tělesové symetrie (sférická / axisymetrická / asymetrická); numerická kontrola zdravého rozumu pro těleso podobné JWST z první verze potvrzující $\tau_{\mathrm{off}} \approx 190\ \text{s}$ vs. integrační okno $600\ \text{s}$ a velikost kroku $0{,}05\ \text{s}$ → 4 řády rezervy; explicitní most k citlivosti momentu síly z gravitačního gradientu pro sestavu dvou těles. Ukotvení: Hughes 2004 §4 + §5 + §8.

---

**Následující kapitola:** [03 — Vzájemný gravitační potenciál a multipolový rozvoj](03-mutual-gravitational-potential.md) *(připravováno)*

**Konec 02-kineticka-energie.md**
