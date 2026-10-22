# Controllo Organic Lab 2.5 → 2.5.1

Data: 22 ottobre 2026. Stato al termine dell’audit: versione preparata localmente. Pubblicazione successivamente **autorizzata dal proprietario il 22 ottobre 2026**.

## Origine e pubblicazione effettiva

Repository: https://github.com/mich090909/organic-lab

Partenza: `640c01c42d4666cc01ecaa2d660796c0af80639b`, ultimo commit di `main` acquisito all'inizio dell'intervento e riconfermato con un nuovo fetch al termine. Quattro file nella radice: `index.html`, `organic-v25-extension.js`, `README.md`, `.nojekyll`. Nessun progetto da ricostruire, servizio backend o dataset esterno. L'HTML contiene il catalogo originale e il percorso guidato; l'estensione aggiunge le voci e i metadati della versione 2.5. L'audit considera il risultato effettivo dell'unione, incluse le precedenti rettifiche legacy.

Il sito pubblico è stato scaricato e confrontato byte per byte con questo commit:

| File | SHA-256 del sorgente e della copia pubblicata |
|---|---|
| index.html | `ab056dbd7cf5516fdd8eb387263421c5bd8e954630e11d86ecebf5046666d12d` |
| organic-v25-extension.js | `0ca72f39444973d48164da37475faebade7d882d16caed4a52bdca9436f71422` |

GitHub conferma la pubblicazione Pages del medesimo commit da `main`: workflow dinamico `pages build and deployment`, [run 39777945598](https://github.com/mich090909/organic-lab/actions/runs/39777945598), concluso con successo. I dati restano nella medesima architettura; nessuna riscrittura dell'app.

## Censimento integrale

**912 voci di nomenclatura + 142 quesiti di isomeria = 1054 voci controllate.** Le molecole riutilizzate nei quesiti non sono altre molecole da sommare al catalogo.

| Famiglia | Molecole |
|---|---:|
| Alcani | 80 |
| Alcheni e alchini | 262 |
| Ciclici e aromatici | 10 |
| Alcoli, fenoli ed eteri | 101 |
| Aldeidi e chetoni | 44 |
| Acidi carbossilici ed esteri | 49 |
| Ammine, ammidi e nitrili | 74 |
| Polifunzionali | 94 |
| Alogenoderivati | 198 |
| **Totale** | **912** |

| Tipo di quesito | Numero |
|---|---:|
| R/S: un centro | 23 |
| R/S: due centri | 17 |
| R/S: tre centri | 10 |
| **Totale assegnazione R/S** | **50** |
| E/Z | 8 |
| Enantiomeri / diastereoisomeri / identici | 54 |
| Isomeria costituzionale / non isomeri | 24 |
| Meso / chirale | 6 |
| **Totale isomeria** | **142** |

Obiettivi numerici già raggiunti nella 2.5: 50 esercizi R/S e oltre 600 molecole. Le 912 voci hanno grafi/configurazioni memorizzate distinti, senza duplicati canonici; alcune rappresentano intenzionalmente una costituzione senza configurazione specificata. La sola presenza di 912 voci non ne dimostrava la correttezza. Sono inoltre presenti 114 Fischer (227 centri complessivi) e 128 alias accettati dal correttore.

## Errori accertati e correzioni

| Problema | Entità e intervento |
|---|---|
| Fischer incoerenti con lo stereoisomero della formula scheletrica | **102 su 114**. Rigenerate tutte le 114 proiezioni dal grafo numerato e dalle priorità complete; ricontrollate leggendo i testi e le coordinate dell'SVG finale. |
| Priorità CIP errate nelle spiegazioni Fischer | **119 centri**, distribuiti nelle stesse 102 proiezioni. Ricalcolate usando il motore CIP moderno, con atomi duplicati per i legami multipli; rimosse le liste basate sull'orientamento/estensione della catena. |
| Numerazione dei nomi | **30 nomi**: 19 con insieme di locanti non minimo, altri 11 risolti dal criterio del sostituente citato per primo quando gli insiemi coincidono. Aggiornati nome, guida, descrittori e quesiti collegati. |
| Prefissi mancanti o sbagliati in Impara | **58 guide originarie**, voci 855–912. Rigenerati dal nome e verificati sul grafo. In totale cambiano 88 campi prefisso, includendo quelli aggiornati per le nuove numerazioni. |
| Conteggio dei centri stereogenici | **151 molecole** contavano soltanto i centri con configurazione già assegnata. Ora il filtro considera tutti i centri effettivi, separatamente dal filtro “stereochimica definita”. |
| Assenza dell'avviso di configurazione non specificata | **131 molecole**. Ora la proprietà è calcolata dai potenziali elementi stereogenici del grafo. |
| Disegni ambigui per stereochimica non specificata | **2** strutture, but-2-ene e pent-2-ene, mostravano una geometria E pur senza specificare E/Z nel nome. Ora il doppio legame è incrociato. |
| Separazione impropria nel nome | **18** nomi avevano un trattino fra l'ultimo prefisso e la radice della catena. Rimossa la separazione; 48 nomi cambiano complessivamente contando anche le 30 rinumerazioni. |
| Locanti e quantità | Uniformato il percorso: posizioni numeriche in ordine crescente, mantenendo le ripetizioni, esclusi i locanti stereochimici e N. Cambiano 140 campi locanti finali; non sono 140 errori chimici, ma una convenzione resa coerente. |
| Precisione delle spiegazioni | **77 testi** precisati: locanti degli alcoli terminali, differenza fra ammine/ammidi non sostituite e N-sostituite, e due diciture locali. Sono stati esaminati tutti i 206 testi distinti originari del catalogo e i modelli dei quesiti. |
| Spiegazione scheletrica | Usava preferenzialmente i metadati Fischer anche quando era mostrata la scheletrica. Ora usa i dati della rappresentazione selezionata e controlla la coerenza del risultato geometrico. |
| Filtro meso | Contava i caratteri `@`, facendo valere `@@` come due centri. Ora usa il conteggio chimico. Corretto anche il conteggio per le coppie. |
| Verifica a tempo | Il codice della 2.5 impostava `endsAt:null` e non avviava il timer. Ripristinati 90 secondi per molecola e l'accesso alla modalità già presente nel codice. |

Le categorie si sovrappongono: **non sommare le righe per ottenere un numero di molecole errate**. La correzione non cambia nessuno dei 912 identificativi o degli stereoisomeri canonici memorizzati. Cinque risposte R/S cambiano i locanti/ordine dopo rinumerazione della catena, senza cambiare lo stereoisomero sottostante.

### Il caso allegato

Nella Fischer originale della voce 733: CH3 in alto; C2 con H a sinistra/F a destra; C3 con Br a sinistra/H a destra; CH2 e CH3 in basso. La decodifica indipendente dà **(2S,3S)**.

Sul C3: **Br > verso C2 [F,C,H] > verso C4 [C,H,H] > H**. Il numero del carbonio o la lunghezza del ramo non stabiliscono la priorità. Con H a destra, davanti all'osservatore, la lettura va invertita.

La formula scheletrica e il nome memorizzati per la stessa voce rappresentavano invece (2S,3R). È stato conservato questo stereoisomero e corretta la Fischer, scambiando Br e H sul C3. Il test di regressione continua a leggere l'originale come (2S,3S): la correzione non si limita a cambiare una lettera della risposta.

## Cosa verifica l'automazione

Il file `inventory.csv` registra tutte le 1054 voci, senza campionamento. `summary.json` riporta conteggi ed eventuali fallimenti.

1. **912 strutture:** parsing e sanitizzazione RDKit, valenze, connettività, formula molecolare, unicità canonica, elementi stereogenici definiti/non definiti.
2. **912 nomi:** ricostruzione separata del grafo a partire dal nome italiano; confronto con SMILES e con il grafo ricostruito dal disegno. Il controllo dei locanti enumera le catene acicliche ammissibili e confronta funzione principale, insaturazioni e sostituenti, compresi i casi di parità alfabetica presenti. I 21 casi contenenti anelli sono elencati in una checklist manuale esplicita e verificati anche per connettività.
3. **912 SVG scheletrici:** lettura dei tracciati visibili per ricostruire atomi, legami, molteplicità, cunei pieni, tratteggi e geometria E/Z. Gli indici degli atomi dell'SVG non sono assunti uguali a quelli degli SMILES canonici. La configurazione è quindi ricavata dalla rappresentazione effettiva, senza usare le lettere della risposta come input geometrico.
4. **114 Fischer:** ricostruzione da etichette e coordinate visibili; orizzontali davanti, verticali dietro. Confronto del grafo e dello stereoisomero con la scheletrica e verifica delle priorità mostrate.
5. **142 soluzioni:** R/S ed E/Z dalla geometria; coppie confrontate per formula, scheletro e funzioni; enantiomeri/identici/diastereoisomeri tramite grafi stereochimici e immagine speculare; meso tramite sovrapponibilità alla propria immagine speculare.
6. **128 alias e 912 guide:** confronto strutturale degli alias sistematici e riferimenti espliciti per i nomi comuni; lunghezza della struttura principale, prefissi, locanti, stereodescrittori e alternative del percorso guidato.

L'oracolo stereochimico è **RDKit 2026.09.1 modern CIPLabeler**, distinto dalle risposte memorizzate nella 2.5 e dal codice dell'interfaccia. Non viene usato il solo `_CIPRank` del vecchio algoritmo. Il codice aggiunto include una routine di riparazione dei dati e un decodificatore separato delle geometrie. Entrambi si appoggiano allo stesso motore chimico: il ricontrollo è indipendente dalle risposte e dal renderer dell'app, **non costituisce un confronto fra due motori cheminformatici diversi**.

## Riferimenti manuali e regressioni

Undici test riproducibili: Fischer allegata, scambio dei due sostituenti di C3, priorità dei due rami carboniosi, rotazioni Fischer 180°/90°, carbonile con duplicazione dell'ossigeno, anello rispetto a catena, prima differenza oltre il primo strato, enantiomeri da SVG, conservazione dei grafi/identificativi, nomi non supportati, precisione delle spiegazioni su N-sostituzione/locanti e checklist finita dei 21 composti con anelli (alcuni aspetti sono raggruppati nello stesso test).

Sono riferimenti controllati sul significato chimico delle strutture e sulle regole, non controlli umani indipendenti di ogni pixel dei 1054 esercizi. La verifica completa delle voci è automatica; l'ispezione visiva manuale riguarda le schermate salvate e i casi di riferimento.

Fonti tecniche: [RDKit, stereochimica](https://www.rdkit.org/docs/RDKit_Book.html#stereochemistry), [RDKit modern CIPLabeler](https://www.rdkit.org/docs/source/rdkit.Chem.rdCIPLabeler.html), [IUPAC P-14, locanti e ordine dei prefissi](https://iupac.qmul.ac.uk/BlueBook/P1.html).

## Interfaccia, filtri e conservazione dei progressi

Il test Playwright usa un browser isolato con uno storico preesistente: conferma che il primo record e la chiave `organic-lab-liceo-v1` restano invariati. Nessun accesso allo storico reale dell'utente. I vecchi risultati non vengono ricalcolati retroattivamente.

Verificate tutte le opzioni del percorso guidato per le 912 molecole, tutte le soluzioni dei 142 quesiti, tutti i nomi e alias nel correttore e le combinazioni dei filtri. Verificati eventi reali dei moduli R/S corretti/errati, cinque passaggi di Impara, prova a tempo, report, esportazione CSV, selezione vuota e comportamento con estensione mancante. Il catalogo incompleto produce un avviso e conserva lo storico. Voci future senza stato di controllo sono escluse dalla selezione.

Su viewport smartphone di 390 px sono state ispezionate le schermate per le voci 733, 797, 836 e 370; verificato anche il layout desktop a 1366 px. Le coppie vengono disposte verticalmente sul telefono; i pulsanti di navigazione vanno a capo e tutti i disegni possono essere ingranditi in una finestra scorrevole.

Il collaudo comprende **3262 rendering**, **2736 insiemi di alternative guidate**, correzione di tutti i **142 quesiti** e **17.160 combinazioni di filtri**, registrati in `ui-summary.json`. Queste sono prove Chromium in emulazione di viewport, **non collaudi fisici su ogni modello di telefono o su Safari iOS**.

## Esito e limiti residui

**Nessun fallimento nei controlli previsti sulle 1054 voci; nessuna voce chimicamente dubbia rimasta nell'archivio corrente e nessuna esclusione necessaria.** Non è una dichiarazione di perfezione universale della nomenclatura.

Il parser dei nomi è limitato al vocabolario italiano presente e agli alias espliciti; non certifica ogni possibile nome preferito IUPAC o sinonimo valido. Il lettore SVG è per i formati effettivi di Organic Lab, non per immagini o disegni arbitrari. Le priorità di nuovi sistemi isotopici, pseudoasimmetrici, non tetraedrici o con geometrie fuori da questi formati richiedono un'estensione e nuovi riferimenti; non devono essere aggiunti con una risposta arbitraria. Il test sui cicli esercita il motore CIP, mentre l'archivio non contiene quesiti R/S ciclici.

La versione è stata preparata e verificata prima dell’approvazione. Il proprietario ha poi autorizzato l’aggiornamento del repository e del sito, includendo codice, dataset, test e rapporto. La cronologia del ramo `main` e i run GitHub Pages documentano la pubblicazione successiva. Le condivisioni non sono state modificate.
