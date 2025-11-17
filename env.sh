#!/bin/bash

# 创建新的 conda 环境
conda create -n flash_dock python=3.11 -y

# 激活环境
conda activate flash_dock

# 先安装 torch 及其相关组件（匹配版本）
echo "正在安装 PyTorch 2.5.1..."

# 尝试清华源（国内更快）
pip install torch==2.5.1 torchvision torchaudio --index-url https://mirrors.tuna.tsinghua.edu.cn/pytorch/whl/cu121

# 验证 torch 版本
TORCH_VERSION=$(python -c "import torch; print(torch.__version__)" 2>/dev/null | cut -d+ -f1)
if [ "$TORCH_VERSION" != "2.5.1" ]; then
    echo "清华源版本不匹配 (安装了 $TORCH_VERSION)，切换到官方源..."
    pip uninstall -y torch torchvision torchaudio
    pip install torch==2.5.1 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121
fi

# 最终验证
python -c "import torch; print(f'✓ PyTorch {torch.__version__} 安装成功')" || {
    echo "❌ PyTorch 安装失败，请检查错误信息"
    exit 1
}

# 卸载可能冲突的旧包
pip uninstall -y xformers datasets transformers 2>/dev/null || true

# 安装 unicore（需要 torch，所以放在 torch 之后）
pip install --no-build-isolation 'unicore @ git+ssh://git@github.com/dptech-corp/Uni-Core.git@ace6fae1c8479a9751f2bb1e1d6e4047427bc134'

# 安装其他依赖
pip install -r requirements.txt

# 重新安装兼容的版本来解决冲突
pip install --upgrade transformers datasets

# 如果需要 xformers，安装兼容版本（可选）
# pip install xformers --no-deps || echo "xformers installation skipped"
