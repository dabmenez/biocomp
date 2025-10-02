#!/bin/bash

# Script simples para configurar ambiente
echo "🧬 Configurando ambiente..."

# Verificar Python
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 não encontrado"
    exit 1
fi

# Criar ambiente virtual
python3 -m venv .venv
source .venv/bin/activate

# Instalar dependências
pip install -r requirements.txt

echo "✅ Ambiente configurado!"
echo "Para ativar: source .venv/bin/activate"