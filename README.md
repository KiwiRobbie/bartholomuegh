# Bartholomuegh
*But you can call me barty!*

![barty](https://github.com/KiwiRobbie/bartholomuegh/assets/46206090/08b70037-2ddb-4096-9f9c-14d6c5568827)

A blackhole rendering project written in rust. 

## Compiling
### Local
```bash
cargo run --release
```

### Web

First, install install `wasm-bindgen` and run the following:
```bash
cargo build --release --target wasm32-unknown-unknown
wasm-bindgen --out-name barty \
  --out-dir web \
  --target web target/wasm32-unknown-unknown/release/bartholomuegh.wasm
```
Then add the following to a html file:
```html
<script type="module">
  import init from './barty.js'
  init()
</script>
```
