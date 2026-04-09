# Structural Dynamics and Aeroelasticity — Assignment 1

| | |
|---|---|
| **Course** | Structural Dynamics and Aeroelasticity |
| **Professors** | Giuseppe Quaranta |
| **Institution** | Politecnico di Milano — Laurea Magistrale Aeronautical Engineering |
| **Date** | April 1, 2026 |
| **Author** | Camile Sam Adam |

---

## Repository Structure

```
.
├── code/
│   ├── src/
│   └── results/
├── Config/
│   ├── background_avion léger1.png
│   ├── background_ulm1.png
│   ├── background_ULM2.png
│   ├── backgroundTemplate.png
│   ├── LogoCN_Q.png
│   ├── Macros.sty
│   ├── Package.sty
│   └── title_page.tex
├── LICENSE-CODE
├── LICENSE-REPORT
├── main.tex
├── README.md
└── report/
    ├── img
    ├── Q1.tex
    ├── Q2.tex
    ├── Q3.tex
    ├── Q4.tex
    ├── Q5.tex
    └── Q6.tex
```

## Build — Report

```bash
latexmk -pdf -interaction=nonstopmode main.tex
```

## Build — Images And Comments

```bash
latexmk -pdf -interaction=nonstopmode main_comments.tex
```

## Build — Numerical Code

```bash
latexmk -pdf -interaction=nonstopmode main_code.tex
```

## License

- **Code** (`code/`): [MIT](LICENSE-CODE)
- **Report** (`report/`): [CC BY-NC 4.0](LICENSE-REPORT)
