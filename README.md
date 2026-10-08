# Organic Lab 2.4 · Quinta liceo scientifico

Web app autonoma per esercitarsi sulla nomenclatura organica e sull’isomeria con 370 molecole e 100 quesiti di isomeria.

## Modalità

- **Impara:** cinque passaggi controllati (struttura principale, locanti, prefissi, stereochimica E/Z e R/S, nome IUPAC finale). Si può ritentare e chiedere un suggerimento; dopo due errori è possibile vedere la soluzione di un passaggio. Per R/S ed E/Z si assegna una configurazione a ogni centro/doppio legame.
- **Allenati:** nome libero e correzione immediata, con suggerimenti facoltativi.
- **Isomeria:** riconoscimento e confronto di isomeri, enantiomeri, diastereoisomeri e meso.

## Pubblicazione

Pubblicare `index.html`, `README.md` e `.nojekyll` nella radice del ramo `main`, poi attivare **Settings → Pages → Deploy from a branch**, `main` e `/ (root)`.

## Privacy e dati

I progressi sono conservati soltanto nel browser dello studente (localStorage: `organic-lab-liceo-v1`); non servono account. Le risposte e le soluzioni sono presenti nel codice client, perciò l’app è destinata all’allenamento, non a verifiche protette. La v2.4 mantiene la chiave dei progressi precedenti e le 370 strutture con relativi nomi.

*Nota:* la valutazione automatica del nome finale usa le denominazioni registrate nel catalogo; non riconosce necessariamente ogni variante IUPAC corretta. Le scelte guidate sono basate sui nomi del catalogo e nei casi più complessi costituiscono una semplificazione didattica della procedura IUPAC.
