# ITU-R S.1503 Notes (itur-s1503-notes)

Notas técnicas e tutoriais (LaTeX + Python) para reproduzir e explicar procedimentos da Recomendação ITU-R S.1503.

## Pré-requisitos

- LaTeX com **LuaLaTeX** (ex.: TeX Live)
- Python 3.10+ (recomendado) e dependências do projeto
- (Opcional) Jupyter para notebooks

## Estrutura

- `notes/001-alpha-angle/`
  - `paper/pt/main.tex` - texto da Nota 001 (Português)
  - `src/alpha_angle.py` - implementação em Python
  - `src/alpha_angle_demo.ipynb` - notebook de demonstração

> Convenção: `notes/NNN-<topic>/` contém `paper/` (LaTeX) e `src/` (código).

## Como compilar a Nota 001

```bash
cd notes/001-alpha-angle/paper/pt
lualatex main.tex
lualatex main.tex
