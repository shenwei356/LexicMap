---
title: Installation
weight: 20
---

LexicMap can be installed via [conda](#conda), [executable binary files](#binary-files),
or [compiling from the source](#compile-from-the-source).

Besides, it supports [shell completion](#shell-completion), which could help accelerate typing.

## Conda

[Install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html), then run

    conda install -c bioconda lexicmap

## Binary files

{{< tabs "uniqueid" >}}
{{< tab "Linux" >}}

1.  Download the binary file.

    |OS     |Arch      |File, 中国镜像                                                                                                                                                                                                               |
    |:------|:---------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
    |Linux  |**64-bit**|[**lexicmap_linux_amd64.tar.gz**](https://github.com/shenwei356/LexicMap/releases/download/v0.4.0/lexicmap_linux_amd64.tar.gz), [中国镜像](http://app.shenwei.me/data/lexicmap/lexicmap_linux_amd64.tar.gz)                  |
    |Linux  |arm64     |[**lexicmap_linux_arm64.tar.gz**](https://github.com/shenwei356/LexicMap/releases/download/v0.4.0/lexicmap_linux_arm64.tar.gz), [中国镜像](http://app.shenwei.me/data/lexicmap/lexicmap_linux_arm64.tar.gz)                  |

2. Decompress it

        tar -zxvf lexicmap_linux_amd64.tar.gz

3. If you have the root privilege, simply copy it to `/usr/local/bin`:

        sudo cp lexicmap /usr/local/bin/

   Or copy to anywhere in the environment variable `PATH`:

        mkdir -p $HOME/bin/; cp lexicmap $HOME/bin/


{{< /tab >}}

{{< tab "MacOS" >}}

1.  Download the binary file.

    |OS     |Arch      |File, 中国镜像                                                                                                                                                                                                               |
    |:------|:---------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
    |macOS  |**64-bit**|[**lexicmap_darwin_amd64.tar.gz**](https://github.com/shenwei356/LexicMap/releases/download/v0.4.0/lexicmap_darwin_amd64.tar.gz), [中国镜像](http://app.shenwei.me/data/lexicmap/lexicmap_darwin_amd64.tar.gz)               |
    |macOS  |arm64     |[**lexicmap_darwin_arm64.tar.gz**](https://github.com/shenwei356/LexicMap/releases/download/v0.4.0/lexicmap_darwin_arm64.tar.gz), [中国镜像](http://app.shenwei.me/data/lexicmap/lexicmap_darwin_arm64.tar.gz)               |


{{< /tab >}}


{{< tab "Windows" >}}

1. Download the binary file.


    |OS     |Arch      |File, 中国镜像                                                                                                                                                                                                               |
    |:------|:---------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
    |Windows|**64-bit**|[**lexicmap_windows_amd64.exe.tar.gz**](https://github.com/shenwei356/LexicMap/releases/download/v0.4.0/lexicmap_windows_amd64.exe.tar.gz), [中国镜像](http://app.shenwei.me/data/lexicmap/lexicmap_windows_amd64.exe.tar.gz)|


2. Decompress it.

2. Copy `lexicmap.exe` to `C:\WINDOWS\system32`.

{{< /tab >}}

{{< tab "Others" >}}

- Please [open an issue](https://github.com/shenwei356/LexicMap/issues) to request binaries for other platforms.
- Or [compiling from the source](#compile-from-the-source).

{{< /tab>}}


{{< /tabs >}}



## Compile from the source


1. [Install go](https://go.dev/doc/install).

        wget https://go.dev/dl/go1.22.4.linux-amd64.tar.gz

        tar -zxf go1.22.4.linux-amd64.tar.gz -C $HOME/

        # or
        #   echo "export PATH=$PATH:$HOME/go/bin" >> ~/.bashrc
        #   source ~/.bashrc
        export PATH=$PATH:$HOME/go/bin

2. Compile LexicMap.

        # ------------- the latest stable version -------------

        go get -v -u github.com/shenwei356/LexicMap/lexicmap

        # The executable binary file is located in:
        #   ~/go/bin/lexicmap
        # You can also move it to anywhere in the $PATH
        mkdir -p $HOME/bin
        cp ~/go/bin/lexicmap $HOME/bin/


        # --------------- the development version --------------

        git clone https://github.com/shenwei356/LexicMap
        cd LexicMap/lexicmap/
        go build

        # The executable binary file is located in:
        #   ./lexicmap
        # You can also move it to anywhere in the $PATH
        mkdir -p $HOME/bin
        cp ./lexicmap $HOME/bin/


## Shell-completion

Supported shell: bash|zsh|fish|powershell

Bash:

    # generate completion shell
    lexicmap autocompletion --shell bash

    # configure if never did.
    # install bash-completion if the "complete" command is not found.
    echo "for bcfile in ~/.bash_completion.d/* ; do source \$bcfile; done" >> ~/.bash_completion
    echo "source ~/.bash_completion" >> ~/.bashrc

Zsh:

    # generate completion shell
    lexicmap autocompletion --shell zsh --file ~/.zfunc/_kmcp

    # configure if never did
    echo 'fpath=( ~/.zfunc "${fpath[@]}" )' >> ~/.zshrc
    echo "autoload -U compinit; compinit" >> ~/.zshrc

fish:

    lexicmap autocompletion --shell fish --file ~/.config/fish/completions/lexicmap.fish
