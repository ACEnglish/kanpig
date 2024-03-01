Building for linux

docker run -v `pwd`:/data --entrypoint /bin/bash -it rust
apt-get update
rustup add target x86_64-unknown-linux-musl
cd /data
cargo build --release --target x86_64-unknown-linux-musl
