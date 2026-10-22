# Organic Lab 2.5 — Palestra di nomenclatura e stereochimica

Web app per la classe quinta del liceo scientifico, utilizzabile dal browser e da smartphone. Non richiede login e non raccoglie dati degli studenti.

**Contenuti verificati (versione 2.5):** 912 strutture stereochimicamente distinte, 142 quesiti di isomeria, tra cui esattamente 50 esercizi per assegnare la configurazione assoluta R/S. I quesiti R/S comprendono 1, 2 o 3 centri stereogenici e sono separati da quelli su enantiomeri, diastereoisomeri e composti meso. Sono presenti 114 proiezioni di Fischer, generate e controllate rispetto alle rispettive formule scheletriche.

## Modalità didattiche

- **Impara**: cinque tappe (struttura principale, locanti, sostituenti, stereochimica, nome IUPAC completo), suggerimenti e correzione a ogni tappa.
- **Allenati**: nome a risposta libera, feedback immediato e spiegazioni.
- **Isomeria**: 50 esercizi R/S con risposte centro per centro, quesiti E/Z, riconoscimento dei rapporti fra stereoisomeri e isomeria costituzionale.

Filtri per famiglia, idrocarburo, difficoltà, numero di centri stereogenici e raffigurazione (scheletrica, Fischer, casuale). Quando selezioni Fischer l'app mostra solo strutture con una Fischer validata; negli altri tipi di esercizi resta il disegno scheletrico. Una struttura non si ripete nella stessa sessione.

## Proiezioni di Fischer

I legami orizzontali puntano verso l'osservatore, i verticali vanno dietro il piano. Solo la rotazione di 180° nel piano mantiene la configurazione; la rotazione di 90° non è equivalente. La correzione mostra le priorità CIP e spiega l'inversione dell'esito quando il gruppo di priorità 4 è davanti.

## Pubblicazione

I seguenti file devono essere pubblicati insieme nella **radice** del ramo `main`: `index.html`, `organic-v25-extension.js`, `README.md`, `.nojekyll`. Attivare GitHub Pages con **Settings → Pages → Deploy from a branch → main → / (root)**.

Indirizzo previsto: https://mich090909.github.io/organic-lab/

**Non eliminare `organic-v25-extension.js`:** contiene la parte ampliata del catalogo ed è necessario per i conteggi completi.

## Privacy e conservazione dei progressi

La chiave `localStorage` resta **`organic-lab-liceo-v1`**, compatibile con le versioni precedenti, e non viene cancellata dal nuovo codice. Lo storico precedente resta sullo stesso browser e dominio. Dati e soluzioni sono incorporati nei file del sito, accessibili pubblicamente: Organic Lab è uno strumento di allenamento, non una verifica protetta. Su dispositivi o origini web differenti i progressi non vengono sincronizzati.

## Accuratezza

I disegni sono prodotti da RDKit; struttura, configurazioni E/Z e R/S e duplicati sono sottoposti a controlli automatici. Sono stati inoltre corretti 100 descrittori stereochimici errati preesistenti nella versione 2.4, mantenendo gli stessi identificativi e le stesse strutture. La nomenclatura segue nomi italiani progettati per le famiglie incluse, e il motore di correzione riconosce i nomi inseriti nel catalogo e le alternative esplicite, **non necessariamente ogni sinonimo IUPAC valido**.
