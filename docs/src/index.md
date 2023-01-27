```@meta
CurrentModule = Main.OurModuleName
DocTestSetup = quote
    println("hi!")
    using Main.OurModuleName
end
```

```@autodocs
Modules = [Main.OurModuleName]
```