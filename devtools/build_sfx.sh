#!/bin/bash
#set -euo pipefail

LIB_DIR=${OPENMM_HOME}/lib
PLUGINS_DIR=${LIB_DIR}/plugins
LIB_FILE=${LIB_DIR}/libOpenMM.so.7.7

function prepare_build() {
    tool_name=$1
    origin=$2
    echo "tool_name: ${tool_name}"
    echo "origin: $origin"

    mkdir ${tool_name}
    cp ${origin}/${tool_name} ${tool_name}
    cp $LIB_FILE ${tool_name}
    cp -r ${PLUGINS_DIR} ${tool_name}/

    cat << EOF > ${tool_name}/run.sh
#!/bin/sh
set -eo pipefail

extracted_dir=\$PWD

echo "Extracted dir: \$extracted_dir"
echo "USER_PWD: \$USER_PWD"

export OPENMM_PLUGIN_DIR=\${extracted_dir}/plugins
export LD_LIBRARY_PATH="\${extracted_dir}":\$LD_LIBRARY_PATH
cd \$USER_PWD
\${extracted_dir}/narupa-cli \$@
EOF
    chmod +x ${tool_name}/run.sh

    ls -ltra ${tool_name}

    makeself $tool_name ${tool_name}.run "$tool_name" ./run.sh
}

target_dir=$(readlink -f target/release)
echo "Target dir is $target_dir"

cargo build --release

build_dir=$(mktemp -d)
pushd ${build_dir}
prepare_build narupa-cli $target_dir
popd
cp ${build_dir}/narupa-cli.run .
cp -r ${build_dir} .

