# Flash Dock 安装指南

## 重要说明

**PyTorch 版本**: 本项目使用 **PyTorch 2.5.1**（当前最新稳定版）。

> ⚠️ 注意：之前的 requirements.txt 中错误地标注了 `torch==2.9.1`，但该版本不存在。已更正为 `2.5.1`。

## 依赖冲突解决方案

### 问题描述

在服务器环境中可能遇到以下依赖冲突：
- `torchaudio`/`torchvision` 版本与 torch 不兼容
- `xformers` 要求特定的 torch 版本
- `transformers` 与新版 `huggingface-hub`/`tokenizers` 冲突
- `datasets` 与新版 `fsspec` 冲突

### 方案 1：全新安装（推荐）

如果是全新环境，直接运行：

```bash
bash env.sh
```

### 方案 2：在已有环境中安装

如果服务器上已有一些预装包，按以下步骤操作：

```bash
# 1. 激活环境
conda activate flash_dock

# 2. 先卸载冲突的包
pip uninstall -y torch torchvision torchaudio xformers transformers datasets tokenizers huggingface-hub

# 3. 重新安装 torch 2.5.1 及配套版本
pip install torch==2.5.1 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# 4. 安装 unicore
pip install --no-build-isolation 'unicore @ git+ssh://git@github.com/dptech-corp/Uni-Core.git@ace6fae1c8479a9751f2bb1e1d6e4047427bc134'

# 5. 安装其他依赖
pip install -r requirements.txt
```

### 方案 3：快速修复当前环境

使用提供的修复脚本：

```bash
bash fix_torch.sh
```

该脚本会自动：
- 卸载错误版本的 torch
- 安装正确版本 (2.5.1)
- 验证安装结果

### 方案 4：忽略依赖冲突（临时方案）

如果只是警告而程序能正常运行，可以继续使用：

```bash
pip install -r requirements.txt --no-deps
# 然后手动安装确实需要的依赖
```

## 验证安装

安装完成后，运行以下命令验证：

```bash
python -c "import torch; import unicore; print(f'Torch: {torch.__version__}'); print('Unicore imported successfully')"
```

## 运行应用

```bash
streamlit run main.py
```

## 常见问题

### Q: 国内网络下载慢
已在 `env.sh` 中配置清华大学镜像源。如果需要切换其他镜像：

**清华源（推荐）：**
```bash
pip config set global.index-url https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple
# PyTorch: --index-url https://mirrors.tuna.tsinghua.edu.cn/pytorch/whl/cu121
```

**阿里云源：**
```bash
pip config set global.index-url https://mirrors.aliyun.com/pypi/simple/
# PyTorch: --index-url https://mirrors.aliyun.com/pytorch-wheels/cu121/
```

**中科大源：**
```bash
pip config set global.index-url https://pypi.mirrors.ustc.edu.cn/simple/
```

**豆瓣源：**
```bash
pip config set global.index-url https://pypi.douban.com/simple/
```

**恢复官方源：**
```bash
pip config unset global.index-url
```

### Q: CUDA 版本不匹配
如果你的 CUDA 版本不是 12.1，需要修改 torch 安装源：

**官方源：**
- CUDA 11.8: `--index-url https://download.pytorch.org/whl/cu118`
- CUDA 12.4: `--index-url https://download.pytorch.org/whl/cu124`
- CPU only: `--index-url https://download.pytorch.org/whl/cpu`

**清华源（国内）：**
- CUDA 11.8: `--index-url https://mirrors.tuna.tsinghua.edu.cn/pytorch/whl/cu118`
- CUDA 12.1: `--index-url https://mirrors.tuna.tsinghua.edu.cn/pytorch/whl/cu121`
- CPU only: `--index-url https://mirrors.tuna.tsinghua.edu.cn/pytorch/whl/cpu`

### Q: SSH 密钥问题
如果无法访问 Uni-Core 仓库，确保：
1. SSH 密钥已添加到 GitHub
2. 或者使用 HTTPS：修改 URL 为 `https://github.com/dptech-corp/Uni-Core.git`

### Q: 内存不足
如果安装过程中内存不足，可以：
```bash
pip install --no-cache-dir -r requirements.txt
```

