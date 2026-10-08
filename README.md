# Organic Lab · Versione studenti 2.3

Palestra di nomenclatura organica per la quinta del liceo scientifico: **370 molecole** e **100 esercizi di isomeria**, comprese configurazioni **E/Z**, **R/S**, molecole con più centri stereogenici e casi misti.

## Avvio

Aprire `index.html` nel browser oppure pubblicare questa cartella tramite **GitHub Pages**.

L'app è un singolo HTML con tutti i disegni SVG incorporati: non richiede compilazione, librerie esterne, connessioni API o account. Offline: scaricare e aprire `index.html` sul dispositivo.

## Funzioni

- `Impara`: suggerimenti e spiegazioni passo per passo.
- `Allenati`: scrittura del nome completo e correzione immediata.
- `Isomeria`: classificazione, E/Z e R/S, enantiomeri, diastereoisomeri e composti meso.
- Filtri per famiglia, sottofamiglia idrocarburica, difficoltà e stereochimica, inclusi centri R/S multipli.
- Statistiche e ripasso degli errori, esportazione CSV e azzeramento dei progressi.

**Non è previsto alcun pannello docente né un sistema di verifica o votazione.** Il catalogo viene aggiornato nel codice e ripubblicato dal proprietario del repository.

## Dati personali e salvataggio

I progressi sono salvati nel `localStorage` del browser con chiave `organic-lab-liceo-v1` e non sono trasmessi a server. Non viene richiesto il nome o una registrazione. I progressi non si sincronizzano tra browser o dispositivi. Non si trasferiscono automaticamente dalla copia `file://` alla versione `https://` di GitHub Pages.

La versione studenti ignora eventuali esclusioni o modifiche al catalogo salvate dal vecchio pannello docente; conserva invece lo storico degli allenamenti sullo stesso sito/origine.

L'intero catalogo con nomi e soluzioni è incorporato nell'HTML ed è quindi leggibile nel sorgente della pagina. Progettato solo per allenamento.

La correzione considera il nome previsto e le alternative elencate, **non tutte le possibili denominazioni sistematiche valide**. Le strutture e gli stereodescrittori sono stati ereditati dalla versione 2.2 e non sono stati modificati in questa revisione dell'interfaccia.

## Pubblicazione GitHub Pages

1. Creare un repository chiamato `organic-lab` nel proprio account GitHub personale.
2. Caricare i file `index.html`, `README.md` e `.nojekyll` nella cartella principale del ramo `main`.
3. Aprire **Settings → Pages → Build and deployment → Deploy from a branch**, selezionare `main` e `/ (root)`, poi salvare.
4. Il sito diventerà raggiungibile all'indirizzo `https://NOME-UTENTE.github.io/organic-lab/`.

**La pubblicazione è pubblica**: non caricare dati degli studenti, credenziali o materiali privati. Il pacchetto qui fornito contiene soltanto il codice dell'app e la documentazione. Nessun file viene pubblicato automaticamente.
