// We need to forces Rust's libpthread linker.
// See https://github.com/extendr/rextendr/issues/77
void R_init_howfar_extendr(void *dll);

// Standard R package entrypoint. R calls this on dyn.load() of the package.
// We forward to the extendr-generated initializer which registers all
// #[extendr]-annotated Rust functions.
void R_init_howfar(void *dll) {
    R_init_howfar_extendr(dll);
}
