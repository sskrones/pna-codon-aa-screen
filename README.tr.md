# PNA kodonu × Bennu amino asidi taraması

*[English version → README.md](README.md)*

Dört kanonik olmayan, prebiyotik açıdan makul bazdan kurulmuş bir peptit nükleik asit (PNA) genetik
sisteminde amino asitler kodonlara **doğrudan stereokimyasal tanıma** ile atanmış olabilir mi?
Bu repo, **64 üç-bazlı PNA kodonunun** tamamını, **Bennu asteroit örneklerinde ve karbonlu göktaşlarında
bildirilen 15 amino aside** (D-enantiyomerler ve protein dışı amino asitler dahil) karşı tarıyor.
Toplam **960 çiftin** her biri AutoDock Vina ile dock edildi, Amber ile implicit solventte simüle edildi ve
MM-GBSA ile puanlandı.

| Kod | Baz | Amber rezidüsü |
|---|---|---|
| **D** | 2,6-diaminopürin | `DAP` |
| **M** | 6-metilurasil | `MTU` |
| **H** | hipoksantin | `HYX` |
| **U** | urasil | `URA` |

Amino asitler: Gly, L-Ala, D-Ala, β-Ala, GABA, Aib (α-aminoizobütirik asit), L-izovalin,
L-Val, D-Val, L-Asp, D-Asp, L-Glu, D-Glu, D-Gln, L-Met.

## Ana sonuç (negatif)

**Bu model kodona özgü bir tanıma olduğuna dair kanıt vermiyor.**

1. **Docking skorları toplamsal ve boyut kaynaklı.** Vina afinitesindeki varyansın %92'si bir amino asit
   terimi (%49) ile bir kodon teriminin (%43) toplamıyla açıklanıyor.
   - Amino asit terimi ligand boyutunu izliyor (ağır atom sayısıyla Spearman ρ = −0,94).
   - Kodon terimi 2,6-diaminopürin sayısını izliyor (ρ = −0,80).
   - Özgüllüğü taşıması gereken kodon × amino asit terimi yalnızca %8 (SD 0,09 kcal/mol). Bu, Vina'nın
     doğruluk sınırının çok altında.
2. **Kompleksler MD'de kararlı değil.** Bütün MM-GBSA ortalamaları sıfırın 2 SD'si içinde
   (−2,5 ile +0,45 kcal/mol arası; %61'i pozitif). Bağımsız OpenMM tekrarlarında ligand onlarca pikosaniye
   içinde trimerden uzaklaşıyor. Yani ΔG büyük ölçüde ligandın ne kadar süre oyalandığını ölçüyor.
   En üst sıradaki çiftler (UMH–Iva −2,48, MDU–Met −2,14 kcal/mol) dörder tohumla tekrarlandığında
   −0,64 ile +0,11 kcal/mol arasında kaldı.
3. **Kiralite tercihi yok.** L/D çiftleri her kodonda benzer dock ediliyor. Yalnızca Glu/D-Glu'da bir fark
   var (medyan 0,05 kcal/mol, p = 0,005) ve bu fark yöntemin doğruluğunun çok altında. MM-GBSA'da L/D
   farkı yok.

**Yorum.** Bu modelde serbest amino asitler ile tek iplikli DMHU-PNA üçlüleri özgül kompleksler oluşturmuyor.
Böyle bir genetik sistem var olsaydı, kodon–amino asit atamasının doğrudan üçlü bağlanmasından çok
adaptör benzeri ara moleküllere ihtiyaç duyması daha olası. Bu sonuç *aşağıda listelenen sınırlar içinde*
geçerli; "tespit edilmedi" ile "imkânsız" aynı şey değil.

| | |
|---|---|
| ![](figures/fig2_variance_decomposition.png) | ![](figures/fig3_aa_size_effect.png) |
| ![](figures/fig5_mmgbsa_ranked.png) | ![](figures/fig7_validation_dissociation.png) |

Tam rapor: [docs/rapor_tr.md](docs/rapor_tr.md) · [docs/report_en.md](docs/report_en.md)

## Klasör yapısı

```
data/
  receptors/   {CODON}.prmtop/.inpcrd (tleap)   {codon}_min.pdbqt (minimize, docking reseptörü)
  ligands/     {aa}.prmtop/.inpcrd/.pdbqt        GAFF2 + AM1-BCC, nötr form
  mdin/        sander girdileri: min_complex, md_equil, md_prod, sp_energy
results/
  mmgbsa_results_md.csv     ham tarama çıktısı, 960 satır
  tables/                   matrisler, özetler, varyans ayrıştırma, kiralite testleri, summary_stats.json
docs/                       tam rapor: rapor_tr.md (Türkçe), report_en.md (English)
figures/                    fig1–fig8 (PNG)
scripts/                    01–07, sırasıyla pipeline (aşağıda)
validation/                 seçili çiftlerin bağımsız OpenMM tekrarları + çizim scripti
notebooks/                  orijinal Colab notebook çıktısı (kaynak kaydı)
```

## Pipeline

| Adım | Script | Ne yapıyor |
|---|---|---|
| 1 | `01_parametrize_pna_monomers.sh` | Kepli PNA monomerleri → GAFF2 tipleri, AM1-BCC yükleri (antechamber, parmchk2) |
| 2 | `02_build_residue_libs.leap` | Kepleri kaldır, head/tail tanımla → zincirlenebilir tleap rezidüleri |
| 3 | `03_build_codon_trimers.py` | 64 trimeri kur, GB minimizasyonu (sander, igb=5) |
| 4 | `04_parametrize_amino_acids.py` | 15 amino asit → GAFF2/AM1-BCC, prmtop/inpcrd |
| 5 | `05_make_pdbqt.py` | GAFF2→SYBYL tip dönüşümü, AM1-BCC yüklü PDBQT, HD/A atom tipleri |
| 6 | `06_dock_md_mmgbsa.py` | Vina docking → ParmEd kompleksi → min / 100 ps ısıtma / 150 ps MD → MM-GBSA (50 kare) |
| 7 | `07_analyze.py` | Tablolar, istatistikler, figür 1–6 |
| — | `validation/rerun_pair_openmm.py` | Aynı protokol OpenMM'de (OBC2, mbondi); ligand–trimer mesafesini izler |

1–5. adımlar zaten yapıldı; çıktıları `data/` içinde. 1. ve 4. adımın ham girdileri (monomer `.cjson`
ve amino asit `.sdf` dosyaları) `data/raw/` klasörüne konacak; bu klasör henüz repoda yok.

### Yeniden üretme

```bash
conda env create -f environment.yml && conda activate pna-screen
python scripts/06_dock_md_mmgbsa.py            # tam tarama (CPU'da saatler sürer, kesilirse kaldığı yerden devam eder)
python scripts/06_dock_md_mmgbsa.py --codons DDD DMH --aas asp gly --out results/test.csv   # alt küme
python scripts/07_analyze.py                    # analiz ve figürler (saniyeler)
python validation/rerun_pair_openmm.py DDD asp 1 > validation/results/seed_DDD_asp_1.json
python validation/plot_validation.py
```

## Bilinen sorunlar ve sınırlar

Bunlar negatif sonucu açıklamıyor; doğrulama tekrarları tutarlı bir kurulumla yapıldı ve aynı davranışı
gösteriyor. Yine de çalışmayı ilerletecek biri bunları düzeltmeli.

- **Reseptör geometrisi uyuşmazlığı.** Docking minimize trimerle (`*_min.pdbqt`) yapılıyor, kompleks ise
  minimize edilmemiş tleap koordinatlarıyla (`{CODON}.inpcrd`) kuruluyor. İkisi arasındaki RMSD 0,4–1,1 Å.
- **Açık trimer uçları.** Rezidü kütüphanesi hazırlanırken kepler kaldırılıyor ama yerine uç grup
  eklenmiyor. Bu yüzden ilk omurga N'si ve son karbonil C'si ikişer bağlı kalıyor.
- **Nötr amino asitler.** Tüm ligandlar net yük 0 (COOH/NH₂) ile kuruldu. pH 7 civarında sudaki gerçek
  formlar zwitteriyonlar ve Asp⁻/Glu⁻ anyonları.
- **GB yarıçapları.** Topolojiler `mbondi` kullanıyor; `igb=5` için önerilen set `mbondi2`.
- **Kısıtlamasız ve kısa MD.** Zayıf bir komplekste 150 ps'lik serbest MD ligandın ayrılmasına izin
  veriyor. Entropisiz tek-trajektori MM-GBSA bu durumda afiniteyi değil, ligandın kalış süresini
  raporluyor.
- **Modelin kapsamı.** Tek iplikli trimerler, implicit solvent, çift başına tek dock pozu, metal iyonu ya da
  mineral yüzey yok, standart dışı PNA rezidüleri için genel bir kuvvet alanı (GAFF2).

## Önerilen sonraki adımlar

1. Yukarıdaki sorunları düzelt: tutarlı reseptör koordinatları, kepli uçlar, zwitteriyonlar, mbondi2.
2. Serbest MD yerine zayıf bağlayıcılara uygun bir serbest enerji yöntemi kullan: kısıtlı MD ile umbrella
   sampling ya da alkimyasal ΔG, veya explicit solventte çift başına çok sayıda kısa replika.
3. Bağlanma ceplerinin oluşabileceği daha büyük yapılı sistemleri test et (PNA dupleksleri, saç tokası
   ya da halka yapıları). Stereokimyasal hipotezin RNA aptamer çalışmaları da bu tür yapılarla yapılıyor.
4. Docking'i yalnızca filtre olarak kullan ve ligand boyutunu açıkça kontrol et (ör. ligand verimliliği).

## Kaynaklar

İngilizce README'ye bakın.

## Yazar

Pelinsu Eda Şentürk, Biyomühendislik, İstanbul Medeniyet Üniversitesi


## Lisans

Kod: MIT (`LICENSE`). Veri ve figürler: CC BY 4.0.
