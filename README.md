# potencial-inhibitors-of-sPLA2-IIA

## Najważniejsze informacje:

Znajdujący się w repozytorium kod jest częścią pracy magisterskiej pt. "Analiza roli fosfolipazy A2 i jej potencjalnych inhibitorów w rozwoju infekcji spowodowanej wirusem SARS-CoV-2" (ang. 'Analysis of the role of phospholipase A2 and its potential inhibitors in the development of infection caused by SARS-CoV-2 virus')

Badania wskazują na znaczenie ludzkiej fosfolipazy A2 podczas sztormu cytokinowego występującego w ostrym przebiegu COVID-19. Zahamowanie tego enzymu może zniwelować ciężkie objawy i zapobiec śmierci pacjenta.

Sekwencje potencjalnych inhibitorów wybrane do analizy: **FLSYK, HDMNKVLDL, LAIYS, VDIHVWDGV (pip9), LGRVDIHVWDGVYIRGR (pip17) i VDIHVWDGVVDIHVWDGV (pip18)**.

Ligandy reprezentujące fosfolipazę A2:

**1bbc_AF2_rel_m1** – wynik zwinięcia programem AF2 sekwencji znajdującej się  w rekordzie 1bbc w bazie PDB; skrót rel oznacza, że model został poddany relaksacji, natomiast m1 to oznaczenia modelu, który w AF2 miał największą wartość  pLDDT, czyli był najwyżej oceniony przez program

**1bbc_AF2_unrel_m1** – jak wyżej; różnica polega na tym, że struktura nie została poddana relaksacji w polu siłowym Amber

**1bbc_PDB** – struktura, którą możemy wyszukać w bazie PDB wpisując kod 1bbc; została pobrana bezpośrednio z bazy

**prekursorowa_rel_m3** – output AF2  powstały w wyniku zwinięcia sekwencji prekursorowej fosfolipazy A2, poddany relaksacji w polu siłowym; jest to model  oznaczony cyfrą 3, ponieważ to właśnie on znajdował się najwyżej w rankingu AF2

**prekursorowa_unrel_m3** – jak wyżej; struktura niepoddana relaksacji

## Streszczenie

Ludzkie fosfolipazy A2, m. in. sPLA2-IIA, są odpowiedzialne za degradację fosfolipidów błonowych, również tych, które budują cząsteczki patogenów, w związku z czym wspomniane enzymy pełnią kluczową rolę w zwalczaniu infekcji. Wszakże podczas ciężkiego przebiegu COVID-19 u pacjentów zaobserwowano zawyżony poziom tego enzymu, co koreluje z burzą cytokinową odpowiedzialną za ostry stan zapalny. Nadaktywność układu immunologicznego podczas zakażenia SARS-CoV-2 przyczynia się do dysfunkcji narządów, co może prowadzić do śmierci pacjenta. W leczeniu COVID-19 kluczowe jest zatem znalezienie odpowiedniego inhibitora właściwej fosfolipazy, aby nie doprowadzić do nadekspresji mediatorów zapalnych. W niniejszej pracy zostały omówione substancje drobnocząsteczkowe oraz peptydy, które są kandydatami na inhibitory sPLA2-IIA. Dodatkowo, dla niektórych krótkich peptydów przeprowadziłam badania in silico z wykorzystaniem języka Python 3.9.7 oraz platform symulacyjnych AlphaFold2 oraz CABS-dock służącymi do przewidywania struktury białek oraz umożliwiającymi dokowanie potencjalnych krótkich białkowych inhibitorów do enzymów. Spośród sześciu peptydów o sekwencji aminokwasów: FLSYK, HDMNKVLDL, LAIYS, VDIHVWDGV, LGRVDIHVWDGVYIRGR i VDIHVWDGVVDIHVWDGV zostały wytypowane dwa (LAIYS, VDIHVWDGV), które na podstawie wyników obliczeń są najlepszymi kandydatami na inhibitory ludzkiej fosfolipazy sPLA2-IIA. Kryteria wykorzystane do oceny badanych związków to zdolność do wiązania z miejscem aktywnym fosfolipazy, z fragmentem sekwencji przyłączającym jon wapniowy oraz regionem odpowiedzialnym za wzmożoną aktywność katalityczną. Skrypty użyte do obliczeń znajdują się w repozytorium Github: https://github.com/Aniczk/potencial-inhibitors-of-sPLA2-IIA. Przyszłe badania mogą polegać na dalszej funkcjonalizacji wytypowanych oligopeptydowych inhibitorów w celu uzyskania ich większej specyficzności, a dalej chemicznej syntezy jak również eksperymentalnego sprawdzenia ich biologicznej aktywności.

## Abstract

Human phospholipases A2, e.g. sPLA2-IIA, are responsible for the degradation of membrane phospholipids, including those that build pathogen molecules, therefore these enzymes play a key role in fighting infection. After all, during the severe course of COVID-19, elevated levels of this enzyme were observed in patients, which correlates with the cytokine storm responsible for acute inflammation. Overactivity of the immune system during SARS-CoV-2 infection contributes to organ dysfunction, which can lead to the death of the patient. In the treatment of COVID-19, it is therefore crucial to find the right inhibitor of the right phospholipase so as not to overexpress inflammatory mediators. This paper discusses small molecule substances and peptides that are candidates for sPLA2-IIA inhibitors. In addition, for some short peptides, I conducted in silico tests using Python 3.9.7 and simulation platforms AlphaFold 2 and CABS-dock for protein structure prediction and for docking potential short protein inhibitors to enzymes. Out of six peptides with the amino acid sequences: FLSYK, HDMNKVLDL, LAIYS, VDIHVWDGV, LGRVDIHVWDGVYIRGR and VDIHVWDGVVDIHVWDGV, two (LAIYS, VDIHVWDGV) were selected, which, based on the calculation results, are the best candidates for the human sPLA2-IIA phospholipase inhibitors. The criteria used to evaluate the tested compounds are the ability to bind to the active site of the phospholipase, to the fragment of the sequence that attaches the calcium ion or to the site responsible for the increased catalytic activity. The scripts used for the calculations are available in the Github repository: https://github.com/Aniczk/potencial-inhibitors-of-sPLA2-IIA. Future research may be related to further functionalization of selected oligopeptide inhibitors in order to obtain their better specificity, and further, in chemical synthesis as well as experimental verification of their biological activity.

