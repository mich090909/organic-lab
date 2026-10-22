# Organic Lab 2.5.1 — Nomenclatura e stereochimica

Web app statica per il liceo, senza account, con progressi locali. Mantiene il catalogo e le funzioni della versione 2.5; corregge le rappresentazioni Fischer, le spiegazioni CIP, la numerazione di alcuni nomi e i metadati dei filtri.

912 voci di nomenclatura; 142 quesiti di isomeria: 50 R/S (23 con un centro, 17 con due, 10 con tre), 8 E/Z, 54 relazioni stereoisomeriche, 24 isomeria costituzionale, 6 meso/chirali. Sono disponibili 114 proiezioni di Fischer.

Il [rapporto di controllo](audit/REPORT.md) distingue controlli automatici, riferimenti manuali, correzioni e limiti. Il [censimento completo](audit/inventory.csv) contiene una riga per ciascuna delle 1054 voci controllate. I conteggi non implicano che ogni voce sia una sostanza con configurazione completamente specificata: alcune strutture sono intenzionalmente prive di stereodescrittori.

## Utilizzo

- **Impara:** cinque passaggi, dalla struttura principale al nome completo.
- **Allenati:** risposta libera, correzione e spiegazione.
- **Verifica:** prova a tempo, senza soluzioni durante lo svolgimento.
- **Isomeria:** assegnazione R/S ed E/Z, coppie di isomeri e forme meso.

Filtri per famiglia, idrocarburo, difficoltà, stereochimica, numero effettivo di centri stereogenici e rappresentazione. Il numero dei centri include quelli la cui configurazione non è specificata; il filtro della stereochimica definita è distinto. I locanti numerici nel percorso guidato sono proposti in ordine crescente, conservando le ripetizioni: indicano posizioni, non quantità.

I pulsanti **Ingrandisci struttura** permettono di leggere i dettagli anche su smartphone. Le formule possono essere ingrandite e scorse senza modificarne la geometria.

## Progressi

La chiave resta `organic-lab-liceo-v1`. Identificativi e strutture molecolari precedenti sono conservati. Lo storico non viene cancellato o riscritto: le vecchie risposte attese possono quindi mostrare i nomi precedenti alla correzione. Nuove sessioni usano i dati corretti. I progressi restano nel browser e sul dominio usato, senza sincronizzazione fra dispositivi. È disponibile l'esportazione CSV.

## Pubblicazione

File necessari, nella radice di `main`: `index.html`, `organic-v25-extension.js`, `.nojekyll`. Aggiornarli insieme; mantenere anche questo README. GitHub Pages pubblica da `main`, cartella `/`, con il workflow dinamico `pages build and deployment`. Non è necessario ricostruire l'app o installare dipendenze sul sito.

Sito: https://mich090909.github.io/organic-lab/

Pubblicazione della versione 2.5.1 autorizzata dal proprietario il 22 ottobre 2026. Per aggiornamenti successivi, verificare che `main` non sia avanzato rispetto alla base di lavoro e integrare eventuali novità. Non forzare aggiornamenti del ramo.

## Controlli riproducibili

Da un clone con la cronologia del repository e Python 3.12:

```sh
python -m pip install -r tests/requirements.txt
python tests/audit_catalog.py
python -m unittest discover -s tests -p 'test_*.py' -v
npm install --no-save playwright@1.62.1
npx playwright install chromium
node tests/ui.cjs
```

`tests/repair_catalog.py` documenta e riproduce la riparazione dei dati dal commit originale fissato in `tests/catalog.py`; preserva le modifiche separate dell'interfaccia. Non va usato per importare nuovi esercizi o aggiornamenti successivi. Nessuno degli script pubblica o modifica repository remoti.

## Perimetro scientifico

Le priorità sono calcolate durante la preparazione con **RDKit 2026.09.1, modern CIPLabeler**, sul grafo completo; l'app non usa lunghezze dei rami o locanti come scorciatoie. La verifica ricostruisce anche i grafi e la stereochimica dalle linee, dai cunei/tratteggi e dai testi degli SVG effettivi. Include legami multipli, anelli, simmetria e confronto con l'immagine speculare.

Il controllo dei nomi comprende il vocabolario italiano presente nel catalogo e 128 alias espliciti, non ogni possibile sinonimo IUPAC. Non è un motore generale di nomenclatura o di riconoscimento di immagini. Gli esercizi futuri senza stato di controllo sono esclusi dalla selezione; sintassi o geometrie non supportate fanno fallire i test. Le soluzioni sono pubbliche nel codice: usare la modalità Verifica come simulazione didattica.
