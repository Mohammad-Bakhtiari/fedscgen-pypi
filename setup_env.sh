#!/bin/bash

set -e  # Exit on any command failure

# ----------- Ensure conda is fully initialized for future shells ------------
echo "🔧 Ensuring conda is initialized for shell..."
if ! grep -q 'conda initialize' ~/.bashrc; then
    echo "Running 'conda init bash' to initialize conda in your shell..."
    conda init bash
    echo "✅ Conda initialization complete. Restart your shell or run 'source ~/.bashrc' for changes to take effect."
else
    echo "ℹ️  Conda already initialized in your shell."
fi

# ----------- Create and activate main environment ---------------------------
echo "📦 Creating environment from environment.yml..."
export SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True
conda env create -f environment.yml || {
    echo "❌ Environment creation failed."
    exit 1
}

echo "🔄 Activating 'fedscgen' environment..."
# Immediate session activation
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate fedscgen || {
    echo "❌ Failed to activate fedscgen environment."
    exit 1
}

echo "✅ Environment 'fedscgen' successfully created and activated."
echo "🚀 Ready to go. You can now run your experiments."

# ----------- R environment setup --------------------------------------------
echo ""
echo "🧪 Creating R environment from r_eval.yml..."
if conda env list | grep -q "r_eval"; then
    echo "ℹ️  R environment 'r_eval' already exists. Skipping creation."
else
    conda env create -f r_eval.yml || {
        echo "❌ R environment creation failed."
        exit 1
    }
fi

echo "🔄 Activating 'r_eval' environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate r_eval || {
    echo "❌ Failed to activate r_eval environment."
    exit 1
}

echo "📦 Installing R libraries via install_libraries.R..."
if Rscript install_libraries.R; then
    echo "✅ R libraries installed successfully."
else
    echo "❌ Failed to install R libraries."
    exit 1
fi

echo "✅ R environment 'r_eval' is ready."
echo ""
echo "🎉 All environments set up successfully!"
