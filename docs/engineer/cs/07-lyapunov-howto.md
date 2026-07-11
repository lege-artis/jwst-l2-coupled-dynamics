# Inženýrská úroveň — návod: výpočet Ljapunovových spekter

> **Publikum.** Praktik, jenž chce *spouštět* nástroje pro Ljapunovovo spektrum v `lyapunov/` proti trajektorii a odečíst odpověď — nikoli odvozovat algoritmus. Pro odvození z prvních principů (BGGS, vnoření $\mathbb{R}^{14}$, kritérium integrability) viz kanonickou úroveň: `docs/canonical/cs/07-ljapunovovo-spektrum.md`.
>
> **Stav.** v0.1 (pokrývá jednovektorové nástroje fáze 2 + nástroje plného spektra fáze 3). CS překlad `docs/engineer/en/07-lyapunov-howto.md`.
>
> **Stručně.** Chcete znaménko chaosu (je $\lambda_{\max} > 0$?)? Užijte `lyapunov_spectrum(..., n_vectors=1)`. Chcete *plné* 12rozměrné spektrum + hamiltonovské párování? Užijte `lyapunov_spectrum_deflated(...)` (levné) nebo `lyapunov_spectrum_liealg(...)` (nativní pro varietu). Neužívejte plochý ovladač $\mathbb{R}^{14}$ při `n_vectors=14` pro plné spektrum — vnáší dva nefyzikální kontrahující módy (viz §5).

---

## §1 Kdy sáhnout po kterém nástroji

| Chcete… | Volejte | Cena |
|---|---|---|
| Pouze $\lambda_{\max}$ (znaménko chaosu) | `lyapunov.spectrum.lyapunov_spectrum(..., n_vectors=1)` | 1 trajektorie + 1 tečný vektor |
| Plné 12rozměrné spektrum, levně | `lyapunov.cocycle_deflated.lyapunov_spectrum_deflated(...)` | 1 trajektorie + 12 tečných vektorů, $\mathbb{R}^{14}$ konečná diference + deflace |
| Plné 12rozměrné spektrum, nativní pro varietu | `lyapunov.cocycle_liealg.lyapunov_spectrum_liealg(...)` | stejný počet; kocyklus $\mathfrak{so}(3)\times\mathbb{R}^3$ |
| Kalibrovat pipeline vůči známé fyzice | `testbeds.lyapunov_three_scenario` (a/b/c) | tři referenční běhy |
| Reprodukovat brány plného spektra fáze 3 | `testbeds.lyapunov_full_spectrum` | ovladače Tier A/B |

Všechny čtyři ovladače sdílejí týž tvar volání a vracejí týž `dict` (viz §4). Krok integrátoru `step_fn` je normálně `integrators.symplectic_se3_variational.se3_variational_step` — symplektický KDK krok, jenž zachovává moment hybnosti na ULP, takže referenční trajektorie je čistá.

## §2 Rychlý start — $\lambda_{\max}$ pro scénář Džanibekova

```python
import numpy as np
from dynamics import BodyState, pack_state
from integrators.symplectic_se3_variational import se3_variational_step
from lyapunov.spectrum import lyapunov_spectrum
from lyapunov.reference import dzhanibekov_characteristic_frequency

# Asymetrické těleso poblíž sedla prostřední osy, primár o hmotnosti Země nablízku.
I_a = np.diag([100.0, 200.0, 400.0])      # I1 < I2 < I3 (nutné pro Džanibekova)
m_a = 500.0
omega0 = np.array([1e-3, 1e-2, 1e-3])     # poblíž e2 (prostřední osa)
m_b = 5.97e24
I_b = np.eye(3) * (2/5 * m_b * 6.371e6**2)
rho = 3.0e6

body_a = BodyState(np.zeros(3), np.zeros(3), np.array([1.,0,0,0]), omega0)
body_b = BodyState(np.array([rho,0,0]), np.zeros(3), np.array([1.,0,0,0]), np.zeros(3))
state0 = pack_state(body_a, body_b)

result = lyapunov_spectrum(
    state0, se3_variational_step, dt=50.0, n_steps=100_000, n_vectors=1,
    m_a=m_a, I_a_body=I_a, m_b=m_b, I_b_body=I_b,
)

lam_max = result['exponents'][0]
sigma   = dzhanibekov_characteristic_frequency(100., 200., 400., omega0[1])
print(f"lambda_max = {lam_max:.3e} s^-1   (sigma = {sigma:.3e})")
# lambda_max ~ 1.5e-5 > 0  -> slabě chaotický (gravitační gradient narušil integrabilitu)
```

Kladné `lam_max` znamená deterministický chaos; `lam_max ≈ 0` s omezeným pohybem znamená regulární/integrabilní dynamiku. Zde ~1,5e-5 > 0 potvrzuje, že vazba gravitačního gradientu činí pohyb poblíž sedla slabě chaotickým.

## §3 Plné 12rozměrné spektrum

Pro celé spektrum (a kontrolu hamiltonovského párování) užijte 12rozměrný kocyklus. Dva ekvivalentní vstupní body — oba berou **tytéž argumenty jako `lyapunov_spectrum`** a mají výchozí `n_vectors=12`:

```python
from lyapunov.cocycle_deflated import lyapunov_spectrum_deflated   # Tier A
from lyapunov.cocycle_liealg  import lyapunov_spectrum_liealg      # Tier B

res = lyapunov_spectrum_deflated(
    state0, se3_variational_step, dt=50.0, n_steps=4_000,
    m_a=m_a, I_a_body=I_a, m_b=m_b, I_b_body=I_b,
)
print(res['exponents'])        # 12 hodnot, sestupně
print(res['sum_exponents'])    # mělo by být << 0.14 (blízko 0 pro hamiltonovský tok)
```

- **Tier A** (`_deflated`) ponechává konečnodiferenční kocyklus $\mathbb{R}^{14}$, avšak před každým krokem Gramovy–Schmidtovy ortonormalizace promítá každý tečný vektor ortogonálně ke dvěma radiálním kvaternionovým směrům (jež se vyvíjejí s trajektorií). Nejlevnější cesta k čistému 12rozměrnému spektru.
- **Tier B** (`_liealg`) propaguje tečné vektory v tělesové Lieově algebře $\mathfrak{so}(3)\times\mathbb{R}^3$, jež je 12rozměrná *konstrukcí* — směr kvaternionové normy se nikdy neobjeví. Stejný počet vyhodnocení (1 + 2·n_vectors volání integrátoru na krok). Užijte, chcete-li nativní odpověď pro varietu bez následné projekce.

Oba dávají spektrum, jež se sčítá na ≈ 0 a podporuje hamiltonovské párování (§6).

## §4 Čtení výsledného slovníku

Každý ovladač vrací tytéž klíče:

| Klíč | Význam |
|---|---|
| `exponents` | `n_vectors` Ljapunovových exponentů, sestupně (s⁻¹) |
| `sum_exponents` | součet vypočtených exponentů (ovladače plného spektra; očekává se ≈ 0) |
| `tau_qr_used` | produkční kadence reortonormalizace QR (s), automaticky laděná |
| `lambda_max_pilot` | odhad $\lambda_{\max}$ z levného pilotního běhu |
| `n_qr` | počet provedených kroků QR |
| `running_exponents` | `[(t, exponents), …]` snímky — vykreslete je k viditelnosti konvergence |
| `T_total` | celkový simulovaný čas v produkčním běhu (s) |
| `final_state` | 26složkový stav na konci |

Kadence QR je automaticky laděná: krátký pilotní běh odhadne $\lambda_{\max}$, pak je produkční kadence nastavena na ~0,1 e-násobku na epochu (ořezáno). Normálně `tau_qr` nenastavujete ručně — předejte `None` (výchozí) a nechte jej naladit.

## §5 Jediná věc, kterou nedělat

**Nevolejte plochý `lyapunov_spectrum(..., n_vectors=14)` a nečtěte jej jako plné spektrum.** Plochá mapa $\mathbb{R}^{14}$ nese dva nedynamické směry kvaternionové normy (po jednom na těleso), jež normalizace v každém kroku integrátoru silně kontrahuje. Objeví se jako dva velké záporné exponenty (≈ −0,064, −0,071 s⁻¹) bez kladných partnerů, takže spektrum $\mathbb{R}^{14}$ se sčítá na ≈ −0,14 místo 0 a **selhává v hamiltonovském párování**. Jednovektorové $\lambda_{\max}$ je nedotčeno (sleduje nejvíce expandující *dynamický* směr), a proto je `n_vectors=1` v pořádku — avšak pro plné spektrum vždy užijte 12rozměrný deflační nebo Lieovsko-algebraický ovladač. Plná diagnóza: kanonická kap. 07 §7.7 + §7.10; empirický záznam `_audit/phase3-prep/FULL-SPECTRUM-FINDING-2026-06-07.md`.

## §6 Brány — kontrola zdravého rozumu spektra

Pomocné rutiny v `lyapunov.reference` a `lyapunov.pairing` převádějí spektrum na kontroly vyhověl/nevyhověl:

```python
from lyapunov.pairing import spectrum_sum_gate, spurious_modes_absent, full_spectrum_gates
from lyapunov.reference import hamiltonian_pairing_check

spectrum_sum_gate(res['exponents'])        # |sum| < 1e-2  (hamiltonovské zachování objemu)
spurious_modes_absent(res['exponents'])    # žádný exponent < -1e-2 (normové módy pryč)
hamiltonian_pairing_check(res['exponents'])# {'pairing_pass': True, 'n_near_zero': ...}
full_spectrum_gates(res['exponents'])      # vše výše sbaleno
```

Pro jednovektorové brány $\lambda_{\max}$ vůči známé fyzice má `lyapunov.reference` scénářové brány:
`lambda_max_gate_scenario_a` (integrabilní: $\lambda_{\max} < 10^{-3}$), `lambda_max_gate_scenario_b` (librace: $|\lambda_{\max}| < \omega_{\mathrm{lib}}$), `lambda_max_gate_scenario_c` (Džanibekov: $\lambda_{\max} > 0$).

## §7 Úskalí

- **Integrabilní případy konvergují pomalu.** Pro regulární/integrabilní trajektorii klesá konečnočasový $\lambda_{\max}$ k 0 jen jako $O(\log T / T)$, takže uvidíte např. ~2e-4 při $T=20\,000$ s spíše než ~0. Toto je **vlastní kvaziperiodické dynamice, nikoli vnoření** — Lieovsko-algebraický ovladač ukazuje týž pomalý pokles (fáze 3 to potvrdila; kanonická kap. 07 §7.3.2/§7.10). Nenastavujte práh integrabilního případu těsněji než ~1e-3 a posuzujte „nulu“ vůči relevantní dynamické frekvenci (σ, $\omega_{\mathrm{lib}}$), nikoli vůči 0.
- **Vykreslete `running_exponents`** předtím, než uvěříte jedinému konečnému číslu — pokud se neustálilo na plató, integrujte déle.
- **Reálný čas** škáluje s n_vectors. 12vektorový běh přes ~4\,000 kroků je desítky sekund; rozpočtujte odpovídajícím způsobem (užijte `production_run=False` v ovladačích testbedu pro rychlou kontrolu).
- **Chaos scénáře (c) je slabý** ($\lambda_c \sim 1{,}5\times10^{-5}$; časová konstanta divergence $1/\lambda \approx 6{,}6\times10^4$ s) — režim KAM poblíž separatrisy, nikoli silná stochasticita. Malý-ale-kladný je očekávanou, správnou odpovědí.

## §8 Kam dál

- **Proč** je algoritmus takový, jaký je, a kritérium integrability: kanonická `docs/canonical/cs/07-ljapunovovo-spektrum.md` (EN: `docs/canonical/en/07-lyapunov-spectrum.md`).
- **Tři kalibrační scénáře** v plném znění: `testbeds/lyapunov_three_scenario.py`.
- **Ovladače + brány plného spektra fáze 3**: `testbeds/lyapunov_full_spectrum.py`, `tests/test_lyapunov_phase3.py`.

---

**Konec 07-lyapunov-howto.md (inženýrská úroveň v0.1, CS).**
